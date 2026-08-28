#!/usr/bin/env python3
"""Stream the exact PLACO+ family and emit canonical independent shared loci."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


VARIANT_FIELDS = {
    "pair_id", "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2",
    "T_PLACO_PLUS", "P_PLACO_PLUS", "analysis_status",
}
SUMMARY_FIELDS = {
    "pair_id", "input_variant_count", "eligible_variant_count", "z2_excluded_count",
    "VarZ1", "VarZ2", "CorZ", "marginal_p_threshold", "analysis_status",
}
CLUMP_FIELDS = {
    "pair_id", "SNP", "locus_id", "lead_snp", "CHR", "start", "end",
    "r2_to_lead", "clumping_reference_id",
}
LOCUS_FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "external_phenotype_name",
    "locus_id", "CHR", "start", "end", "lead_snp", "lead_A1", "lead_A2",
    "lead_Z_sleep", "lead_Z_external", "lead_P_sleep", "lead_P_external",
    "lead_T_PLACO_PLUS", "lead_P_PLACO_PLUS", "primary_locked_family_threshold",
    "conventional_genomewide_status", "sleep_effect_direction", "external_effect_direction",
    "CorZ", "VarZ1", "VarZ2", "significant_variant_count", "significant_variants",
    "clumping_reference_id", "clumping_parameters", "pleiotropy_class", "claim_limit",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def probability(value: str, field: str, identity: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field} for {identity}") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise SystemExit(f"ERROR: invalid {field} for {identity}")
    return number


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/pleiotropy_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/pleiotropy_manifest.lock.json"))
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--variant-results", type=Path, required=True)
    parser.add_argument("--clump-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/pleiotropy/novel_shared_loci.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/pleiotropy_results.json"))
    args = parser.parse_args()

    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: PLACO+ manifest differs from its pre-result lock")
    pair_ids = [row["pair_id"] for row in manifest]
    if pair_ids != lock.get("pair_ids_in_locked_order"):
        raise SystemExit("ERROR: PLACO+ pair family/order differs from lock")
    manifest_by_pair = {row["pair_id"]: row for row in manifest}
    for row in manifest:
        for path_field, checksum_field in (
            ("merged_genomewide_path", "merged_genomewide_sha256"),
            ("ld_reference_path", "ld_reference_sha256"),
            ("placo_source_path", "placo_source_sha256"),
        ):
            path = Path(row[path_field])
            if not path.is_file() or sha256(path) != row[checksum_field]:
                raise SystemExit(f"ERROR: locked PLACO+ input drifted: {row['pair_id']}/{path_field}")

    summary_rows = read_tsv(args.run_summary)
    if not summary_rows or not SUMMARY_FIELDS.issubset(summary_rows[0]) or [row["pair_id"] for row in summary_rows] != pair_ids:
        raise SystemExit("ERROR: PLACO+ run summary is not the exact locked pair family")
    summary_by_pair = {row["pair_id"]: row for row in summary_rows}
    for row in summary_rows:
        identity = row["pair_id"]
        if row["analysis_status"] != "PLACO_PLUS_COMPLETE":
            raise SystemExit(f"ERROR: PLACO+ pair scan did not complete: {identity}")
        if int(row["eligible_variant_count"]) != int(manifest_by_pair[identity]["eligible_variant_count"]):
            raise SystemExit(f"ERROR: eligible PLACO+ variant count differs from lock: {identity}")
        if float(row["marginal_p_threshold"]) != float(lock["nuisance_marginal_p_threshold"]):
            raise SystemExit(f"ERROR: PLACO+ nuisance threshold differs from lock: {identity}")
        for field in ("VarZ1", "VarZ2"):
            if not math.isfinite(float(row[field])) or float(row[field]) <= 0:
                raise SystemExit(f"ERROR: invalid {field}: {identity}")
        if not math.isfinite(float(row["CorZ"])) or not -1 < float(row["CorZ"]) < 1:
            raise SystemExit(f"ERROR: invalid CorZ: {identity}")

    threshold = float(lock["primary_locked_family_threshold"])
    pair_index = {pair_id: index for index, pair_id in enumerate(pair_ids)}
    counts: Counter[str] = Counter()
    failure_counts: Counter[str] = Counter()
    hits: dict[tuple[str, str], dict[str, str]] = {}
    last_pair_index = -1
    last_coordinate: tuple[int, int, str] | None = None
    with args.variant_results.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not VARIANT_FIELDS.issubset(reader.fieldnames or []):
            raise SystemExit(f"ERROR: PLACO+ variant result lacks fields: {sorted(VARIANT_FIELDS - set(reader.fieldnames or []))}")
        for row in reader:
            pair_id = row["pair_id"]
            if pair_id not in pair_index:
                raise SystemExit(f"ERROR: PLACO+ variant row outside lock: {pair_id}")
            current_pair_index = pair_index[pair_id]
            if current_pair_index < last_pair_index:
                raise SystemExit("ERROR: PLACO+ variant results are not in locked pair order")
            coordinate = (int(row["CHR"]), int(row["BP"]), row["SNP"])
            if not 1 <= coordinate[0] <= 22 or coordinate[1] < 1:
                raise SystemExit(f"ERROR: invalid PLACO+ coordinate: {pair_id}/{row['SNP']}")
            if current_pair_index == last_pair_index and last_coordinate is not None and coordinate <= last_coordinate:
                raise SystemExit(f"ERROR: PLACO+ rows are unsorted or duplicated: {pair_id}/{row['SNP']}")
            if current_pair_index != last_pair_index:
                last_coordinate = None
            last_pair_index, last_coordinate = current_pair_index, coordinate
            counts[pair_id] += 1
            identity = f"{pair_id}/{row['SNP']}"
            for field in ("P1", "P2"):
                probability(row[field], field, identity)
            if row["analysis_status"] == "PLACO_PLUS_COMPLETE":
                placo_p = probability(row["P_PLACO_PLUS"], "P_PLACO_PLUS", identity)
                for field in ("Z1", "Z2", "T_PLACO_PLUS"):
                    if not math.isfinite(float(row[field])):
                        raise SystemExit(f"ERROR: invalid {field} for {identity}")
                if float(row["Z1"]) ** 2 > float(lock["primary_z_squared_maximum"]) or float(row["Z2"]) ** 2 > float(lock["primary_z_squared_maximum"]):
                    raise SystemExit(f"ERROR: primary PLACO+ result retained Z squared above lock: {identity}")
                if placo_p < threshold:
                    hits[(pair_id, row["SNP"])] = row
            elif row["analysis_status"] == "FAILED_NUMERICAL_INTEGRATION":
                failure_counts[pair_id] += 1
            else:
                raise SystemExit(f"ERROR: invalid PLACO+ variant analysis status: {identity}")
    for pair_id in pair_ids:
        if counts[pair_id] != int(summary_by_pair[pair_id]["eligible_variant_count"]):
            raise SystemExit(f"ERROR: PLACO+ result row count differs from summary/lock: {pair_id}")

    clumps = read_tsv(args.clump_map)
    if clumps and not CLUMP_FIELDS.issubset(clumps[0]):
        raise SystemExit(f"ERROR: PLACO+ clump map lacks fields: {sorted(CLUMP_FIELDS - set(clumps[0]))}")
    clump_keys = {(row["pair_id"], row["SNP"]) for row in clumps}
    if len(clump_keys) != len(clumps) or clump_keys != set(hits):
        raise SystemExit("ERROR: PLACO+ clump map is not the exact primary-significant variant family")
    grouped: dict[tuple[str, str], list[tuple[dict[str, str], dict[str, str]]]] = defaultdict(list)
    for clump in clumps:
        pair_id = clump["pair_id"]
        if clump["clumping_reference_id"] != manifest_by_pair[pair_id]["ld_reference_id"]:
            raise SystemExit(f"ERROR: PLACO+ clumping reference differs from lock: {pair_id}")
        if int(clump["CHR"]) != int(hits[(pair_id, clump["SNP"])]["CHR"]):
            raise SystemExit(f"ERROR: PLACO+ clump chromosome mismatch: {pair_id}/{clump['SNP']}")
        grouped[(pair_id, clump["locus_id"])].append((hits[(pair_id, clump["SNP"])], clump))

    loci: list[dict[str, object]] = []
    for (pair_id, locus_id), entries in grouped.items():
        lead_ids = {clump["lead_snp"] for _, clump in entries}
        if len(lead_ids) != 1:
            raise SystemExit(f"ERROR: PLACO+ locus has inconsistent lead SNP: {pair_id}/{locus_id}")
        lead_id = next(iter(lead_ids))
        lead_matches = [result for result, _ in entries if result["SNP"] == lead_id]
        if len(lead_matches) != 1:
            raise SystemExit(f"ERROR: PLACO+ lead SNP absent/duplicated: {pair_id}/{locus_id}")
        lead = lead_matches[0]
        clump_rows = [clump for _, clump in entries]
        bounds = {(row["CHR"], row["start"], row["end"], row["clumping_reference_id"]) for row in clump_rows}
        if len(bounds) != 1:
            raise SystemExit(f"ERROR: PLACO+ locus bounds/reference inconsistent: {pair_id}/{locus_id}")
        chromosome, start, end, reference_id = next(iter(bounds))
        manifest_row = manifest_by_pair[pair_id]
        summary = summary_by_pair[pair_id]
        loci.append({
            "pair_id": pair_id, "sleep_trait": manifest_row["sleep_trait"],
            "extension_trait_id": manifest_row["extension_trait_id"],
            "external_phenotype_name": manifest_row["external_phenotype_name"],
            "locus_id": locus_id, "CHR": chromosome, "start": start, "end": end,
            "lead_snp": lead_id, "lead_A1": lead["A1"], "lead_A2": lead["A2"],
            "lead_Z_sleep": lead["Z1"], "lead_Z_external": lead["Z2"],
            "lead_P_sleep": lead["P1"], "lead_P_external": lead["P2"],
            "lead_T_PLACO_PLUS": lead["T_PLACO_PLUS"], "lead_P_PLACO_PLUS": lead["P_PLACO_PLUS"],
            "primary_locked_family_threshold": threshold,
            "conventional_genomewide_status": "PASS" if float(lead["P_PLACO_PLUS"]) < float(lock["conventional_variant_threshold"]) else "FAIL",
            "sleep_effect_direction": "POSITIVE" if float(lead["Z1"]) > 0 else ("NEGATIVE" if float(lead["Z1"]) < 0 else "ZERO"),
            "external_effect_direction": "POSITIVE" if float(lead["Z2"]) > 0 else ("NEGATIVE" if float(lead["Z2"]) < 0 else "ZERO"),
            "CorZ": summary["CorZ"], "VarZ1": summary["VarZ1"], "VarZ2": summary["VarZ2"],
            "significant_variant_count": len(entries),
            "significant_variants": ";".join(sorted(result["SNP"] for result, _ in entries)),
            "clumping_reference_id": reference_id, "clumping_parameters": manifest_row["clumping_parameters"],
            "pleiotropy_class": "STATISTICAL_PLEIOTROPY_LOCUS",
            "claim_limit": "NOT_SHARED_CAUSAL_VARIANT_OR_CAUSAL_MECHANISM",
        })
    loci.sort(key=lambda row: (pair_index[str(row["pair_id"])], int(row["CHR"]), int(row["start"]), str(row["locus_id"])))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=LOCUS_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(loci)
    provenance = {
        "schema_version": "1.0.0", "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifest_sha256": sha256(args.manifest), "lock_sha256": sha256(args.lock),
        "run_summary_sha256": sha256(args.run_summary), "variant_results_sha256": sha256(args.variant_results),
        "clump_map_sha256": sha256(args.clump_map), "pair_scan_count": len(pair_ids),
        "eligible_variant_counts": dict(counts), "numerical_failure_counts": dict(failure_counts),
        "primary_significant_variant_count": len(hits), "independent_locus_count": len(loci),
        "primary_locked_family_threshold": threshold, "output": str(args.out), "output_sha256": sha256(args.out),
        "warning": "Statistical pleiotropic association does not establish a shared causal variant or distinguish horizontal from vertical pleiotropy.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PLEIOTROPY_RESULTS_OK pairs={len(pair_ids)} significant_variants={len(hits)} loci={len(loci)}")


if __name__ == "__main__":
    main()
