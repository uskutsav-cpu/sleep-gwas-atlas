#!/usr/bin/env python3
"""Audit or materialize the 16 full-resolution GWAS required downstream."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import fine_mapping_contract
import variant_map


READINESS_FIELDS = [
    "trait_id", "route", "source_id", "raw_path", "raw_bytes", "raw_status",
    "dependency_status", "output_path", "output_status", "blocker",
]
EXPECTED_SCRIPT_PATHS = {
    "scripts/101_prepare_dense_harmonization.py",
    "scripts/01_harmonize.py",
    "scripts/variant_map.py",
    "scripts/liftover_chain.py",
    "scripts/fine_mapping_contract.py",
}
EXPECTED_DENSE_TRAITS = ["parkinson", "mdd", "ibd", "crohn", "uc", "stroke", "longevity"]
EXPECTED_LIFTOVER_TRAITS = ["ms", "asthma", "t2d", "cad", "telomere_length", "melanoma"]
EXPECTED_DIRECT_TRAITS = ["ldl", "hdl", "triglycerides"]
HEX64 = set("0123456789abcdef")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=READINESS_FIELDS, delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
        path.chmod(0o644)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        os.replace(temporary, path)
        path.chmod(0o644)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_policy(root: Path, relative: str) -> tuple[Path, dict[str, object]]:
    path = root / relative
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"dense harmonization policy is unreadable: {exc}")
    if (
        policy.get("schema_version") != "sleep-atlas-dense-harmonization.1"
        or policy.get("analysis_panel") != "atlas-v1.0"
        or policy.get("expected_trait_count") != 16
        or policy.get("minimum_free_bytes_before_materialization") != 20 * 1024**3
        or not isinstance(policy.get("script_sha256"), dict)
        or set(policy.get("script_sha256", {})) != EXPECTED_SCRIPT_PATHS
        or policy.get("dense_map_traits") != EXPECTED_DENSE_TRAITS
        or policy.get("liftover_traits") != EXPECTED_LIFTOVER_TRAITS
        or policy.get("direct_hg19_traits") != EXPECTED_DIRECT_TRAITS
        or policy.get("output_directory") != "data/harmonized_mixer_full"
        or policy.get("readiness_path")
        != "results/tables/dense_harmonization_readiness.tsv"
        or policy.get("provenance_path")
        != "results/tables/dense_harmonization.provenance.json"
    ):
        fail("dense harmonization policy differs from the locked production family")
    locked_paths = {
        "analysis_panel_path": "analysis_panel_sha256",
        "source_registry_path": "source_registry_sha256",
        "schema_registry_path": "schema_registry_sha256",
        "dense_map_policy_path": "dense_map_policy_sha256",
        "dense_mapping_plans_path": "dense_mapping_plans_sha256",
        "liftover_plans_path": "liftover_plans_sha256",
    }
    for path_key, hash_key in locked_paths.items():
        candidate = Path(str(policy.get(path_key, "")))
        if candidate.is_absolute() or ".." in candidate.parts:
            fail(f"unsafe dense harmonization dependency path: {path_key}")
        full = root / candidate
        if not full.is_file() or sha256(full) != policy.get(hash_key):
            fail(f"dense harmonization dependency differs from policy: {candidate}")
    for relative, expected_hash in policy["script_sha256"].items():
        candidate = Path(relative)
        if (
            candidate.is_absolute() or ".." in candidate.parts
            or not (root / candidate).is_file() or sha256(root / candidate) != expected_hash
        ):
            fail(f"dense harmonization script differs from policy: {relative}")
    return path, policy


def publish_provenance(
    root: Path,
    policy_path: Path,
    policy: dict[str, object],
    readiness: list[dict[str, object]],
) -> None:
    readiness_path = root / str(policy["readiness_path"])
    if not readiness_path.is_file() or len(readiness) != 16 or any(
        row["output_status"] != "READY" for row in readiness
    ):
        fail("dense harmonization provenance requires exactly 16 ready outputs")
    outputs = {}
    for row in readiness:
        trait = str(row["trait_id"])
        data = root / str(row["output_path"])
        qc = data.with_name(f"{trait}.qc.txt")
        if not data.is_file() or not qc.is_file():
            fail(f"dense harmonization output family is incomplete: {trait}")
        outputs[trait] = {
            "data_path": str(data.relative_to(root)),
            "data_bytes": data.stat().st_size,
            "data_sha256": sha256(data),
            "qc_path": str(qc.relative_to(root)),
            "qc_bytes": qc.stat().st_size,
            "qc_sha256": sha256(qc),
            "route": row["route"],
            "source_id": row["source_id"],
        }
    payload = {
        "schema_version": "sleep-atlas-dense-harmonization-provenance.1",
        "analysis_panel": policy["analysis_panel"],
        "policy_path": str(policy_path.relative_to(root)),
        "policy_sha256": sha256(policy_path),
        "readiness_path": str(readiness_path.relative_to(root)),
        "readiness_sha256": sha256(readiness_path),
        "script_sha256": policy["script_sha256"],
        "trait_count": 16,
        "outputs": outputs,
        "claim_limit": policy["claim_limit"],
    }
    provenance_path = root / str(policy["provenance_path"])
    if provenance_path.exists():
        try:
            observed = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            fail(f"dense harmonization provenance is unreadable: {exc}")
        if observed != payload:
            fail("immutable dense harmonization provenance differs from current outputs")
        return
    atomic_json(provenance_path, payload)
    print(f"DENSE_HARMONIZATION_SEALED out={provenance_path.relative_to(root)}")


def qc_route_matches(
    root: Path,
    trait: str,
    raw: Path,
    qc: Path,
    policy: dict[str, object],
    map_by_trait: dict[str, dict[str, str]],
    liftover_by_trait: dict[str, dict[str, str]],
) -> bool:
    """Require the QC ledger to describe this exact selected-source route."""
    value = fine_mapping_contract.qc_value
    infile_sha256 = value(qc, "infile_sha256")
    common = (
        value(qc, "trait") == trait
        and value(qc, "infile") == str(raw.absolute())
        and len(infile_sha256) == 64
        and not set(infile_sha256).difference(HEX64)
        and value(qc, "prefilter_strategy") == "not supplied"
    )
    if not common:
        return False
    if trait in policy["dense_map_traits"]:
        plan = map_by_trait[trait]
        expected_map = (root / plan["map_path"]).absolute()
        return (
            value(qc, "variant_map") == str(expected_map)
            and value(qc, "variant_map_strategy") == plan["strategy"]
            and value(qc, "variant_map_scope")
            == "GENOME_WIDE_IMPUTED_VARIANT_IDENTITY"
            and value(qc, "liftover_chain") == "not supplied"
        )
    if trait in policy["liftover_traits"]:
        plan = liftover_by_trait[trait]
        expected_chain = (root / plan["chain_path"]).absolute()
        return (
            value(qc, "variant_map") == "not supplied"
            and value(qc, "variant_map_strategy") == "not supplied"
            and value(qc, "liftover_chain") == str(expected_chain)
            and value(qc, "liftover_chain_sha256") == plan["chain_sha256"]
            and value(qc, "source_build_argument") == "hg38"
        )
    return (
        value(qc, "variant_map") == "not supplied"
        and value(qc, "variant_map_strategy") == "not supplied"
        and value(qc, "liftover_chain") == "not supplied"
        and value(qc, "source_build_argument") == "hg19"
    )


def dense_map_status(root: Path, policy: dict[str, object]) -> tuple[bool, dict[str, object]]:
    map_policy = json.loads(
        (root / str(policy["dense_map_policy_path"])).read_text(encoding="utf-8")
    )
    output = root / str(map_policy["map_path"])
    provenance_path = root / str(map_policy["map_provenance_path"])
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        validated = variant_map.validate_provenance(output)
    except (OSError, UnicodeError, json.JSONDecodeError, variant_map.VariantMapError):
        return False, {}
    ready = (
        provenance == validated
        and provenance.get("schema_version") == map_policy["map_schema_version"]
        and provenance.get("map_scope") == map_policy["map_scope"]
        and provenance.get("map_bytes") == output.stat().st_size
        and provenance.get("map_sha256") == sha256(output)
        and provenance.get("mapped_rows", 0) >= map_policy["minimum_mapped_rows"]
    )
    return ready, provenance if ready else {}


def build_plan(root: Path, policy: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, object]]:
    panel = read_tsv(root / str(policy["analysis_panel_path"]))
    panel_by_trait = {row["trait_id"]: row for row in panel}
    if len(panel) != 45 or len(panel_by_trait) != 45:
        fail("analysis panel no longer contains exactly 45 unique traits")
    source_rows = read_tsv(root / str(policy["source_registry_path"]))
    source_by_id = {row["source_id"]: row for row in source_rows}
    map_rows = read_tsv(root / str(policy["dense_mapping_plans_path"]))
    map_by_trait = {row["trait_id"]: row for row in map_rows}
    liftover_rows = read_tsv(root / str(policy["liftover_plans_path"]))
    liftover_by_trait = {row["trait_id"]: row for row in liftover_rows}
    map_policy = json.loads(
        (root / str(policy["dense_map_policy_path"])).read_text(encoding="utf-8")
    )
    dense_traits = list(policy["dense_map_traits"])
    lifted_traits = list(policy["liftover_traits"])
    direct_traits = list(policy["direct_hg19_traits"])
    all_traits = dense_traits + lifted_traits + direct_traits
    if (
        len(all_traits) != policy["expected_trait_count"]
        or len(set(all_traits)) != len(all_traits)
        or len(source_rows) != len(source_by_id)
        or len(map_rows) != len(map_by_trait)
        or len(liftover_rows) != len(liftover_by_trait)
        or set(map_by_trait) != set(dense_traits)
        or not set(lifted_traits).issubset(liftover_by_trait)
        or not set(all_traits).issubset(panel_by_trait)
    ):
        fail("dense harmonization trait/route family differs from policy")
    for trait in dense_traits:
        plan, panel_row = map_by_trait[trait], panel_by_trait[trait]
        expected_strategy = "BY_RSID_ALLELES" if trait == "mdd" else "BY_COORD_ALLELES"
        if (
            plan["source_id"] != panel_row["source_id"]
            or plan["strategy"] != expected_strategy
            or plan["map_path"] != map_policy["map_path"]
            or plan["map_scope"] != map_policy["map_scope"]
            or plan["output_build"] != "GRCh37/hg19"
            or (trait == "mdd" and panel_row["build"] != "UNRESOLVED")
            or (trait != "mdd" and panel_row["build"] != "hg19")
        ):
            fail(f"dense identity-mapping route differs from policy: {trait}")
    for trait in lifted_traits:
        plan, panel_row = liftover_by_trait[trait], panel_by_trait[trait]
        if (
            plan["source_id"] != panel_row["source_id"]
            or plan["strategy"] != "UCSC_CHAIN_POINT"
            or plan["chain_path"] != "ref/hg38ToHg19.over.chain.gz"
            or plan["source_build"] != "GRCh38/hg38"
            or plan["output_build"] != "GRCh37/hg19"
            or panel_row["build"] != "hg38"
        ):
            fail(f"dense liftover route differs from policy: {trait}")
    for trait in direct_traits:
        if panel_by_trait[trait]["build"] != "hg19":
            fail(f"dense direct-GRCh37 route differs from policy: {trait}")
    map_ready, map_provenance = dense_map_status(root, policy)
    readiness = []
    output_dir = root / str(policy["output_directory"])
    chain_ready: dict[str, bool] = {}
    for trait in lifted_traits:
        lift = liftover_by_trait[trait]
        chain = root / lift["chain_path"]
        chain_ready[trait] = (
            chain.is_file() and chain.stat().st_size == int(lift["chain_bytes"])
            and sha256(chain) == lift["chain_sha256"]
        )
    for trait in all_traits:
        panel_row = panel_by_trait[trait]
        source_id = panel_row["source_id"]
        if source_id not in source_by_id:
            fail(f"selected source is absent from public source registry: {source_id}")
        raw = root / "data/raw" / panel_row["raw_file"]
        data = output_dir / f"{trait}.harmonized.tsv.gz"
        qc = output_dir / f"{trait}.qc.txt"
        if trait in dense_traits:
            route = "GENOME_WIDE_GRCH37_IDENTITY_MAP"
            dependency_ready = map_ready
        elif trait in lifted_traits:
            route = "PINNED_HG38_TO_HG19_POINT_LIFTOVER"
            dependency_ready = chain_ready[trait]
        else:
            route = "DIRECT_GRCH37"
            dependency_ready = True
        output_ready = (
            data.is_file() and data.stat().st_size > 0 and qc.is_file()
            and fine_mapping_contract.qc_is_full_resolution(qc)
            and qc_route_matches(
                root, trait, raw, qc, policy, map_by_trait, liftover_by_trait,
            )
        )
        blockers = []
        if not raw.is_file() or raw.stat().st_size == 0:
            blockers.append("full selected-source raw file is absent")
        if not dependency_ready:
            blockers.append("route dependency is absent or differs from its checksum seal")
        if not output_ready:
            blockers.append("full-resolution canonical output is absent")
        readiness.append({
            "trait_id": trait, "route": route, "source_id": source_id,
            "raw_path": str(raw.relative_to(root)),
            "raw_bytes": raw.stat().st_size if raw.is_file() else 0,
            "raw_status": "READY" if raw.is_file() and raw.stat().st_size else "BLOCKED",
            "dependency_status": "READY" if dependency_ready else "BLOCKED",
            "output_path": str(data.relative_to(root)),
            "output_status": "READY" if output_ready else "BLOCKED",
            "blocker": "; ".join(blockers),
        })
    return readiness, map_by_trait, liftover_by_trait, map_provenance


def materialize_trait(
    root: Path,
    policy: dict[str, object],
    trait: str,
    map_by_trait: dict[str, dict[str, str]],
    liftover_by_trait: dict[str, dict[str, str]],
    map_provenance: dict[str, object],
) -> None:
    panel_by_trait = {
        row["trait_id"]: row
        for row in read_tsv(root / str(policy["analysis_panel_path"]))
    }
    panel = panel_by_trait[trait]
    raw = root / "data/raw" / panel["raw_file"]
    if not raw.is_file() or raw.stat().st_size == 0:
        fail(f"full selected-source raw file is absent for {trait}: {raw}")
    output_dir = root / str(policy["output_directory"])
    output_dir.mkdir(parents=True, exist_ok=True)
    final_data = output_dir / f"{trait}.harmonized.tsv.gz"
    final_qc = output_dir / f"{trait}.qc.txt"
    if final_data.exists() or final_qc.exists():
        fail(f"dense harmonization output already exists for {trait}; no overwrite permitted")
    free = shutil.disk_usage(output_dir).free
    if free < int(policy["minimum_free_bytes_before_materialization"]):
        fail(
            f"dense harmonization requires {policy['minimum_free_bytes_before_materialization']} "
            f"free bytes; observed {free}"
        )
    with tempfile.TemporaryDirectory(prefix=f"dense-{trait}-", dir=output_dir) as temporary:
        temporary_path = Path(temporary)
        command = [
            sys.executable, str(root / "scripts/01_harmonize.py"),
            "--trait", trait, "--config", str(root / str(policy["analysis_panel_path"])),
            "--infile", str(raw), "--outdir", str(temporary_path),
        ]
        if trait in policy["dense_map_traits"]:
            plan = map_by_trait[trait]
            map_path = root / plan["map_path"]
            if not map_provenance:
                fail("genome-wide dense variant map is not sealed")
            command.extend([
                "--variant-map", str(map_path),
                "--variant-map-strategy", plan["strategy"],
                "--expected-variant-map-bytes", str(map_provenance["map_bytes"]),
                "--expected-variant-map-sha256", str(map_provenance["map_sha256"]),
            ])
            if trait != "mdd":
                command.extend(["--source-build", "hg19"])
        elif trait in policy["liftover_traits"]:
            plan = liftover_by_trait[trait]
            command.extend([
                "--source-build", "hg38", "--liftover-chain", str(root / plan["chain_path"]),
                "--expected-liftover-chain-bytes", plan["chain_bytes"],
                "--expected-liftover-chain-sha256", plan["chain_sha256"],
            ])
        else:
            command.extend(["--source-build", "hg19"])
        result = subprocess.run(command, cwd=root, text=True)
        if result.returncode:
            fail(f"dense harmonization failed for {trait}")
        staged_data = temporary_path / final_data.name
        staged_qc = temporary_path / final_qc.name
        if (
            not staged_data.is_file() or staged_data.stat().st_size == 0
            or not staged_qc.is_file()
            or not fine_mapping_contract.qc_is_full_resolution(staged_qc)
        ):
            fail(f"dense harmonization output failed full-resolution validation: {trait}")
        os.replace(staged_data, final_data)
        os.replace(staged_qc, final_qc)
        final_data.chmod(0o644)
        final_qc.chmod(0o644)
    print(f"DENSE_HARMONIZATION_PUBLISHED trait={trait} bytes={final_data.stat().st_size}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/dense_harmonization_policy.json")
    parser.add_argument("--materialize", action="store_true")
    parser.add_argument("--trait", action="append", default=[])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--no-write-readiness", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, policy = load_policy(root, args.policy)
    if args.report_only and (args.materialize or args.trait or args.all):
        fail("--report-only cannot be combined with materialization modes")
    readiness, map_by_trait, liftover_by_trait, map_provenance = build_plan(root, policy)
    known_traits = [row["trait_id"] for row in readiness]
    if args.all and args.trait:
        fail("use either --all or one or more --trait arguments")
    selected = known_traits if args.all else args.trait
    if len(selected) != len(set(selected)) or not set(selected).issubset(known_traits):
        fail("selected dense harmonization traits are duplicated or outside the locked 16")
    if args.materialize and not selected:
        fail("--materialize requires --all or at least one --trait")
    if selected and not args.materialize:
        fail("--trait/--all require --materialize")
    for trait in selected:
        row = next(record for record in readiness if record["trait_id"] == trait)
        prerequisites = row["raw_status"] == "READY" and row["dependency_status"] == "READY"
        if not prerequisites:
            fail(f"dense harmonization prerequisites are blocked for {trait}: {row['blocker']}")
        if row["output_status"] == "READY":
            fail(f"dense harmonization is already complete for {trait}; no overwrite permitted")
        materialize_trait(
            root, policy, trait, map_by_trait, liftover_by_trait, map_provenance,
        )
        readiness, map_by_trait, liftover_by_trait, map_provenance = build_plan(root, policy)
    if not args.no_write_readiness:
        readiness_path = root / str(policy["readiness_path"])
        atomic_tsv(readiness_path, readiness)
    ready = sum(row["output_status"] == "READY" for row in readiness)
    if ready == 16 and not args.report_only and not args.no_write_readiness:
        publish_provenance(root, policy_path, policy, readiness)
    print(f"DENSE_HARMONIZATION_{'READY' if ready == 16 else 'BLOCKED'} ready={ready}/16")
    for row in readiness:
        if row["output_status"] != "READY":
            print(f"BLOCKED {row['trait_id']}: {row['blocker']}")
    return 0 if ready == 16 or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
