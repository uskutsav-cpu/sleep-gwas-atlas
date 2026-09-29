"""Build a checksum-bound descriptive direction table for admitted PLACO candidates.

This is a partial-family, post-selection description of the two harmonized Z
statistics at the same variant. It is not LD signed-r analysis, fine-mapping,
an independent-locus test, or a new hypothesis test.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
PAIRS = ("insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz")
CANDIDATE_PATH = ROOT / "brain6/results/loci/placo_candidate_variants_partial.tsv"
MASTER_PATH = ROOT / "brain6/results/placo/placo_master.tsv"
FAMILY_LOCK_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
PAIR_LOCK_PATH = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
OUT_PATH = ROOT / "brain6/results/loci/placo_candidate_directions_partial.tsv"
PROVENANCE_PATH = ROOT / "brain6/results/loci/placo_candidate_directions_partial.provenance.json"
OUT_FIELDS = [
    "pair_id", "sleep_trait", "brain_disorder", "SNP", "CHR", "BP",
    "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "Z_SLEEP",
    "Z_DISORDER", "direction_relation", "candidate_status", "lead_SNP",
    "r2_to_lead", "locus_id",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def direction_relation(z_sleep: float, z_disorder: float) -> str:
    """Describe signs for the same harmonized effect allele at one variant."""
    if z_sleep == 0 or z_disorder == 0:
        return "ZERO_Z_UNDEFINED"
    return "CONCORDANT" if (z_sleep > 0) == (z_disorder > 0) else "OPPOSING"


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def canonical_float(value: str) -> float:
    result = float(value)
    if not (-float("inf") < result < float("inf")):
        raise ValueError(f"Expected finite numeric value, received {value!r}")
    return result


def build(out_path: Path = OUT_PATH, provenance_path: Path = PROVENANCE_PATH) -> dict[str, Any]:
    candidates = tsv(CANDIDATE_PATH)
    master_rows = tsv(MASTER_PATH)
    master = {row["pair_id"]: row for row in master_rows}
    if len(master) != len(master_rows) or set(master) != set(PAIRS):
        raise ValueError("PLACO master must contain exactly the four admitted partial-family pairs")
    if not candidates:
        raise ValueError("Candidate table is empty; no evidence-backed empty-family result is available")

    candidate_by_pair: dict[str, dict[str, dict[str, str]]] = {pair: {} for pair in PAIRS}
    for row in candidates:
        pair, snp = row["pair_id"], row["SNP"]
        if pair not in candidate_by_pair or snp in candidate_by_pair[pair]:
            raise ValueError(f"Unknown pair or duplicate candidate key: {pair}/{snp}")
        candidate_by_pair[pair][snp] = row

    pair_lock = json.loads(PAIR_LOCK_PATH.read_text(encoding="utf-8"))
    pair_definitions = {
        row["pair_id"]: row for row in pair_lock["pairs"] if row["pair_id"] in PAIRS
    }
    if set(pair_definitions) != set(PAIRS):
        raise ValueError("Reviewed pair lock does not cover the exact partial PLACO set")

    matched: dict[tuple[str, str], dict[str, str]] = {}
    pair_source_hashes: dict[str, str] = {}
    for pair in PAIRS:
        pair_master = master[pair]
        source = ROOT / "brain6/results/placo" / pair / "variants.tsv.gz"
        expected_hash = pair_master["output_sha256"]
        observed_hash = sha256(source)
        if observed_hash != expected_hash or str(source.relative_to(ROOT)) != pair_master["output_path"]:
            raise ValueError(f"PLACO variant artifact differs from the validated master for {pair}")
        pair_source_hashes[str(source.relative_to(ROOT))] = observed_hash
        wanted = candidate_by_pair[pair]
        with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                snp = row["SNP"]
                if snp not in wanted:
                    continue
                key = (pair, snp)
                if key in matched:
                    raise ValueError(f"Duplicate PLACO variant key: {pair}/{snp}")
                candidate = wanted[snp]
                if row["status"] != "TESTED" or row["headline"] != "True":
                    raise ValueError(f"Candidate is not a tested PLACO headline variant: {pair}/{snp}")
                if (row["CHR"], row["BP"]) != (candidate["CHR"], candidate["BP"]):
                    raise ValueError(f"Candidate coordinates changed in PLACO output: {pair}/{snp}")
                if canonical_float(row["P_PLACO"]) != canonical_float(candidate["P_PLACO"]):
                    raise ValueError(f"Candidate P-value differs from PLACO output: {pair}/{snp}")
                matched[key] = row
        if len([key for key in matched if key[0] == pair]) != len(wanted):
            raise ValueError(f"Only {len([key for key in matched if key[0] == pair])}/{len(wanted)} candidates joined for {pair}")

    output_rows = []
    pair_trait = {
        pair: (pair_definitions[pair]["sleep_trait"], pair_definitions[pair]["disease_trait"])
        for pair in PAIRS
    }
    for pair in PAIRS:
        for snp, candidate in sorted(candidate_by_pair[pair].items()):
            row = matched[(pair, snp)]
            z_sleep, z_disorder = canonical_float(row["Z1"]), canonical_float(row["Z2"])
            sleep, disorder = pair_trait[pair]
            output_rows.append({
                "pair_id": pair, "sleep_trait": sleep, "brain_disorder": disorder,
                "SNP": snp, "CHR": row["CHR"], "BP": row["BP"],
                "P_PLACO": row["P_PLACO"], "Q_WITHIN_PAIR": row["Q_WITHIN_PAIR"],
                "Q_BONFERRONI_ACROSS_PAIRS": row["Q_BONFERRONI_ACROSS_PAIRS"],
                "Z_SLEEP": row["Z1"], "Z_DISORDER": row["Z2"],
                "direction_relation": direction_relation(z_sleep, z_disorder),
                "candidate_status": candidate["candidate_status"],
                "lead_SNP": candidate["lead_SNP"], "r2_to_lead": candidate["r2_to_lead"],
                "locus_id": candidate["locus_id"],
            })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    direction_counts: dict[str, dict[str, int]] = {}
    candidate_status_counts: dict[str, dict[str, int]] = {}
    for pair in PAIRS:
        counts = Counter(row["direction_relation"] for row in output_rows if row["pair_id"] == pair)
        direction_counts[pair] = dict(sorted(counts.items()))
        statuses = Counter(row["candidate_status"] for row in output_rows if row["pair_id"] == pair)
        candidate_status_counts[pair] = dict(sorted(statuses.items()))
    aggregate = dict(sorted(Counter(row["direction_relation"] for row in output_rows).items()))
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6-placo-partial-candidate-direction-v1",
        "status": "PASS_PARTIAL_FAMILY_DESCRIPTIVE",
        "interpretation": (
            "Post-selection sign description at the same harmonized variant. Z1 and Z2 are aligned "
            "to trait 1's allele by the checksum-pinned pair harmonizer. Concordant/opposing signs "
            "are not signed LD, independent replication, causal evidence, or a new hypothesis test. "
            "Only four of five selected tracks are present; 24 candidates lack exact LD-reference positions."
        ),
        "scope": {"pair_ids": list(PAIRS), "n_pairs": len(PAIRS), "n_candidates": len(output_rows),
                  "complete_five_track_family": False},
        "direction_counts": {"all": aggregate, "by_pair": direction_counts},
        "candidate_status_counts": candidate_status_counts,
        "methods": {
            "direction_rule": "CONCORDANT iff sign(Z_SLEEP) == sign(Z_DISORDER); otherwise OPPOSING; zero Z is undefined",
            "sign_basis": "Z1/Z2 are the pair join's harmonized Z statistics on trait 1 A1; B effects are aligned to A1",
            "multiplicity": "No new test or correction; source PLACO P, within-pair q, and across-pair q are carried through unchanged",
        },
        "sources": {
            "builder_script": {"path": str(Path(__file__).resolve().relative_to(ROOT)),
                               "sha256": sha256(Path(__file__).resolve())},
            "candidate_table": {"path": str(CANDIDATE_PATH.relative_to(ROOT)), "sha256": sha256(CANDIDATE_PATH)},
            "placo_master": {"path": str(MASTER_PATH.relative_to(ROOT)), "sha256": sha256(MASTER_PATH)},
            "placo_family_lock": {"path": str(FAMILY_LOCK_PATH.relative_to(ROOT)), "sha256": sha256(FAMILY_LOCK_PATH)},
            "pair_lock": {"path": str(PAIR_LOCK_PATH.relative_to(ROOT)), "sha256": sha256(PAIR_LOCK_PATH)},
            "pair_variant_files": pair_source_hashes,
            "pair_harmonizer": {"path": "extensions/brain6/brain6/gwas.py", "sha256": sha256(ROOT / "extensions/brain6/brain6/gwas.py")},
        },
        "output": {"path": str(out_path.relative_to(ROOT)), "sha256": sha256(out_path), "rows": len(output_rows)},
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    provenance["provenance_sha256"] = sha256(provenance_path)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT_PATH)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE_PATH)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
