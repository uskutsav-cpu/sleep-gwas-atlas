#!/usr/bin/env python3
"""Build the fixed, LD-pruned CATlas variant universe and overlap cache."""
from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

from liftover_chain import load_chain
from catlas_policy import reference_policy_sha256


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(root: Path, value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        fail(f"unsafe {label} path: {value}")
    return root / relative


def require_pin(path: Path, size: int, digest: str, label: str) -> None:
    if not path.is_file() or path.stat().st_size != size or sha256(path) != digest:
        fail(f"{label} is absent or differs from its exact pin")


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    temporary.replace(path)


def read_ids(path: Path) -> list[str]:
    values = [line.strip().split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not values or len(values) != len(set(values)):
        fail("PLINK pruning output is empty or contains duplicate variants")
    return values


def load_bim(path: Path, selected: set[str]) -> dict[str, tuple[int, int]]:
    observed: dict[str, tuple[int, int]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 6:
                fail(f"reference BIM has a malformed row at line {line_number}")
            snp = fields[1]
            if snp not in selected:
                continue
            try:
                chromosome, position = int(fields[0]), int(fields[3])
            except ValueError:
                fail(f"reference BIM has an invalid selected coordinate at line {line_number}")
            if not 1 <= chromosome <= 22 or position <= 0 or snp in observed:
                fail(f"reference BIM has an invalid or duplicate selected variant: {snp}")
            observed[snp] = (chromosome, position)
    if set(observed) != selected:
        fail("PLINK pruning output contains variants absent from the reference BIM")
    return observed


def load_maf(path: Path, selected: set[str], minimum: float) -> dict[str, float]:
    observed: dict[str, float] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().split()
        if not {"SNP", "MAF"}.issubset(header):
            fail("PLINK frequency output lacks SNP and MAF")
        for line in handle:
            values = line.split()
            if len(values) != len(header):
                fail("PLINK frequency output contains a malformed row")
            row = dict(zip(header, values))
            snp = row["SNP"]
            if snp not in selected:
                fail(f"PLINK frequency output contains an unpruned variant: {snp}")
            try:
                maf = float(row["MAF"])
            except (TypeError, ValueError):
                fail(f"PLINK frequency output contains invalid MAF for {snp}")
            if not math.isfinite(maf) or not minimum <= maf <= 0.5 or snp in observed:
                fail(f"PLINK frequency output contains out-of-policy or duplicate MAF for {snp}")
            observed[snp] = maf
    if set(observed) != selected:
        fail("PLINK frequency output does not cover the complete pruned universe")
    return observed


def mark_interval(
    positions_by_chromosome: dict[int, list[tuple[int, str]]],
    starts_by_chromosome: dict[int, list[int]], chromosome: str, start: int, end: int,
) -> list[str]:
    if not chromosome.startswith("chr"):
        return []
    try:
        chrom = int(chromosome[3:])
    except ValueError:
        return []
    positions = positions_by_chromosome.get(chrom, [])
    starts = starts_by_chromosome.get(chrom, [])
    index = bisect_left(starts, start + 1)
    hits: list[str] = []
    while index < len(positions) and positions[index][0] <= end:
        hits.append(positions[index][1])
        index += 1
    return hits


def annotate_catlas(
    ccre_path: Path, peaks_path: Path, selected_cells: list[dict[str, str]],
    mapped: dict[str, tuple[int, int]], prefix: str,
) -> tuple[set[str], dict[str, set[str]], dict[str, int]]:
    positions_by_chromosome: dict[int, list[tuple[int, str]]] = {}
    for snp, (chromosome, position) in mapped.items():
        positions_by_chromosome.setdefault(chromosome, []).append((position, snp))
    for values in positions_by_chromosome.values():
        values.sort()
    starts_by_chromosome = {
        chromosome: [position for position, _ in values]
        for chromosome, values in positions_by_chromosome.items()
    }
    adult: set[str] = set()
    with gzip.open(ccre_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"#Chromosome", "Start", "End", "Present in adult tissues"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fail("CATlas cCRE universe header differs from the locked schema")
        for row in reader:
            if row["Present in adult tissues"] != "yes":
                continue
            try:
                start, end = int(row["Start"]), int(row["End"])
            except ValueError:
                fail("CATlas cCRE universe contains a nonnumeric interval")
            adult.update(mark_interval(
                positions_by_chromosome, starts_by_chromosome,
                row["#Chromosome"], start, end,
            ))

    membership: dict[str, set[str]] = {snp: set() for snp in mapped}
    cell_counts: dict[str, int] = {}
    with zipfile.ZipFile(peaks_path) as archive:
        names = {member.filename for member in archive.infolist() if not member.is_dir()}
        for cell in selected_cells:
            member = cell["archive_member"]
            if member not in names:
                fail(f"selected CATlas peak member is absent: {member}")
            cell_id = f"{prefix}::{cell['metadata_cell_type']}"
            hits: set[str] = set()
            with archive.open(member) as raw, gzip.GzipFile(fileobj=raw) as compressed:
                text = io.TextIOWrapper(compressed, encoding="ascii")
                for line in text:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 3:
                        fail(f"selected CATlas peak member contains a malformed row: {member}")
                    try:
                        start, end = int(fields[1]), int(fields[2])
                    except ValueError:
                        fail(f"selected CATlas peak member contains a nonnumeric interval: {member}")
                    hits.update(mark_interval(
                        positions_by_chromosome, starts_by_chromosome,
                        fields[0], start, end,
                    ))
            hits.intersection_update(adult)
            cell_counts[cell_id] = len(hits)
            for snp in hits:
                membership[snp].add(cell_id)
    return adult, membership, cell_counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["catlas_adult_v4"]
    manifest_path = safe_relative(root, spec["component_manifest"], "CATlas manifest")
    if not manifest_path.is_file():
        fail("CATlas manifest is absent")
    require_pin(manifest_path, manifest_path.stat().st_size, spec["component_manifest_sha256"], "CATlas manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    for component in manifest["components"]:
        path = safe_relative(root, component["path"], f"CATlas {component['component_id']}")
        require_pin(path, int(component["bytes"]), component["sha256"], f"CATlas {component['component_id']}")
        components[component["component_id"]] = path

    cache_path = safe_relative(root, spec["variant_cache_path"], "CATlas variant cache")
    provenance_path = safe_relative(root, spec["variant_cache_provenance_path"], "CATlas cache provenance")
    if cache_path.exists() or provenance_path.exists():
        fail("CATlas fixed-universe cache is immutable and already exists")

    reference_prefix = safe_relative(root, spec["analysis_reference_prefix"], "analysis reference")
    reference_provenance_path = reference_prefix.with_suffix(".provenance.json")
    if not reference_provenance_path.is_file():
        fail("analysis reference provenance is absent")
    reference_provenance = json.loads(reference_provenance_path.read_text(encoding="utf-8"))
    reference_paths: dict[str, Path] = {}
    for suffix in (".bed", ".bim", ".fam"):
        path = reference_prefix.with_suffix(suffix)
        pin = reference_provenance.get("artifacts", {}).get(path.name, {})
        require_pin(path, int(pin.get("bytes", -1)), str(pin.get("sha256", "")), f"analysis reference {suffix}")
        reference_paths[suffix] = path

    runtime_path = root / policy["causal_inference"]["component_manifest"]
    if not runtime_path.is_file():
        fail("causal runtime manifest is absent")
    require_pin(
        runtime_path, runtime_path.stat().st_size,
        policy["causal_inference"]["component_manifest_sha256"], "causal runtime manifest",
    )
    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    plink_pin = next(
        row for row in runtime["components"]
        if row["component_id"] == "PLINK_1_9_STABLE_UNIVERSAL_BINARY"
    )
    plink = safe_relative(root, plink_pin["path"], "PLINK binary")
    require_pin(plink, int(plink_pin["bytes"]), plink_pin["sha256"], "PLINK binary")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_directory = Path(tempfile.mkdtemp(prefix=".catlas-reference-", dir=cache_path.parent))
    try:
        prune_prefix = temporary_directory / "prune"
        prune_command = [
            str(plink), "--bfile", str(reference_prefix), "--maf", str(spec["minimum_maf"]),
            "--indep-pairwise", f"{spec['ld_pruning_window_kb']}kb", "1",
            str(spec["ld_pruning_r2"]), "--allow-no-sex", "--out", str(prune_prefix),
        ]
        pruned = subprocess.run(prune_command, cwd=root, capture_output=True, text=True, check=False)
        prune_in = Path(str(prune_prefix) + ".prune.in")
        if pruned.returncode or not prune_in.is_file():
            fail(f"PLINK CATlas LD pruning failed: {(pruned.stdout + pruned.stderr).strip()}")
        pruned_ids = read_ids(prune_in)
        selected = set(pruned_ids)

        frequency_prefix = temporary_directory / "frequency"
        frequency_command = [
            str(plink), "--bfile", str(reference_prefix), "--extract", str(prune_in),
            "--freq", "--allow-no-sex", "--out", str(frequency_prefix),
        ]
        frequency = subprocess.run(frequency_command, cwd=root, capture_output=True, text=True, check=False)
        frequency_path = Path(str(frequency_prefix) + ".frq")
        if frequency.returncode or not frequency_path.is_file():
            fail(f"PLINK CATlas frequency run failed: {(frequency.stdout + frequency.stderr).strip()}")
        coordinates = load_bim(reference_paths[".bim"], selected)
        maf = load_maf(frequency_path, selected, float(spec["minimum_maf"]))

        chain_spec = policy["regulatory_build_harmonization"]
        chain, chain_provenance = load_chain(
            root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"],
        )
        mapped: dict[str, tuple[int, int]] = {}
        liftover_counts: dict[str, int] = {}
        for snp in pruned_ids:
            status, target = chain.map_point(*coordinates[snp])
            liftover_counts[status] = liftover_counts.get(status, 0) + 1
            if status == "mapped" and target is not None:
                mapped[snp] = (target[0], target[1])
        if not mapped:
            fail("no pruned CATlas reference variant uniquely lifts to GRCh38")

        adult, memberships, cell_counts = annotate_catlas(
            components["CCRE_UNIVERSE"], components["CELL_TYPE_RESTRICTED_PEAKS"],
            manifest["selected_cells"], mapped, spec["cell_type_id_prefix"],
        )
        cell_order = [
            f"{spec['cell_type_id_prefix']}::{cell['metadata_cell_type']}"
            for cell in manifest["selected_cells"]
        ]
        rows: list[dict[str, object]] = []
        for snp in sorted(mapped, key=lambda value: (*coordinates[value], value)):
            chromosome38, position38 = mapped[snp]
            rows.append({
                "SNP": snp, "CHR37": coordinates[snp][0], "BP37": coordinates[snp][1],
                "MAF": format(maf[snp], ".12g"), "CHR38": chromosome38, "BP38": position38,
                "IN_ADULT_CCRE": int(snp in adult),
                "CELL_TYPE_IDS": ";".join(cell_id for cell_id in cell_order if cell_id in memberships[snp]) or "NA",
            })
        atomic_gzip_tsv(cache_path, spec["variant_cache_fields"], rows)
        provenance = {
            "schema_version": "sleep-atlas-catlas-reference.2",
            "reference_policy_sha256": reference_policy_sha256(policy),
            "component_manifest_sha256": sha256(manifest_path),
            "component_sha256": {
                identity: sha256(path) for identity, path in sorted(components.items())
            },
            "analysis_reference_provenance_sha256": sha256(reference_provenance_path),
            "analysis_reference_artifacts": {
                path.name: sha256(path) for path in reference_paths.values()
            },
            "plink_version": pruned.stdout.splitlines()[0] if pruned.stdout else "UNKNOWN",
            "prune_command": prune_command,
            "frequency_command": frequency_command,
            "prune_log_sha256": sha256(Path(str(prune_prefix) + ".log")),
            "frequency_log_sha256": sha256(Path(str(frequency_prefix) + ".log")),
            "chain": chain_provenance,
            "counts": {
                "reference_variants": int(reference_provenance["counters"]["allele_matched_variants"]),
                "ld_pruned_variants": len(pruned_ids),
                "unique_grch38_lifted_variants": len(mapped),
                "adult_ccre_background_variants": len(adult),
                "selected_cell_types": len(cell_order),
            },
            "liftover_status_counts": liftover_counts,
            "cell_background_variant_counts": cell_counts,
            "cache_path": str(cache_path.relative_to(root)),
            "cache_bytes": cache_path.stat().st_size,
            "cache_sha256": sha256(cache_path),
            "claim_limit": spec["claim_limit"],
            "script_sha256": sha256(Path(__file__).resolve()),
        }
        atomic_json(provenance_path, provenance)
    finally:
        shutil.rmtree(temporary_directory, ignore_errors=True)
    print(
        f"CATLAS_REFERENCE_OK variants={len(rows)} adult_ccre={len(adult)} "
        f"cells={len(cell_order)} out={cache_path.relative_to(root)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
