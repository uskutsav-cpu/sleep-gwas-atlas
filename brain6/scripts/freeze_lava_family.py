#!/usr/bin/env python3
"""Freeze the Brain6 local-rg test family before any Brain6 LAVA results exist."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
from brain6.io import atomic_text, json_hash, sha256, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="v1", choices=["v1", "v2"])
    parser.add_argument("--overlap-dir", type=Path,
                        help="Pair-specific LDSC intercept-derived LAVA covariance family")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.out is None:
        args.out = ROOT / f"brain6/config/lava_family_{args.version}.json"
    if args.overlap_dir is not None:
        args.overlap_dir = args.overlap_dir.resolve()
    selection = ROOT / "brain6/config/deep_tracks_v1.tsv"
    selection_sha_path = selection.with_suffix(selection.suffix + ".sha256")
    pair_lock = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
    loci = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
    checksum_path = args.out.with_suffix(args.out.suffix + ".sha256")
    if args.out.exists() or checksum_path.exists():
        raise SystemExit("Refusing to replace an existing local-rg family lock")

    lock = json.loads(pair_lock.read_text(encoding="utf-8"))
    if lock.get("synthetic", False):
        raise SystemExit("Refusing to define a production family from a synthetic pair lock")
    lock_body = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if json_hash(lock_body) != lock.get("lock_sha256"):
        raise SystemExit("Reviewed pair-lock checksum is invalid")

    selection_expected = selection_sha_path.read_text(encoding="utf-8").split()[0]
    selection_actual = sha256(selection)
    if selection_expected != selection_actual:
        raise SystemExit("Frozen deep-track selection checksum is invalid")
    with selection.open(encoding="utf-8", newline="") as stream:
        selected = list(csv.DictReader(stream, delimiter="\t"))
    pairs = [f"{row['primary_sleep_trait']}__{row['disorder']}"
             for row in selected if row.get("eligibility") == "POST_ATLAS_PRIMARY"]
    if len(pairs) != 5 or set(pairs) != {pair["pair_id"] for pair in lock["pairs"]}:
        raise SystemExit("Local-rg pair family differs from the reviewed five-pair lock")
    if not loci.is_file():
        raise SystemExit(f"Pinned LAVA locus definition is missing: {loci}")
    locus_sha = sha256(loci)
    if locus_sha != "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882":
        raise SystemExit("Pinned LAVA locus-definition checksum mismatch")
    with loci.open(encoding="utf-8") as stream:
        locus_rows = sum(1 for _ in stream)
    if locus_rows != 2496:
        raise SystemExit(f"Expected header plus 2,495 loci, observed {locus_rows} lines")

    overlap = None
    if args.overlap_dir is not None:
        overlap_manifest = args.overlap_dir / f"lava_overlap_family_{args.version}.tsv"
        overlap_rows = list(csv.DictReader(overlap_manifest.open(encoding="utf-8", newline=""),
                                           delimiter="\t"))
        if (len(overlap_rows) != len(pairs)
                or {row["pair_id"] for row in overlap_rows} != set(pairs)
                or any((row["trait1"], row["trait2"]) != next(
                    (pair["sleep_trait"], pair["disease_trait"])
                    for pair in lock["pairs"] if pair["pair_id"] == row["pair_id"])
                       for row in overlap_rows)
                or any(row["pair_lock_sha256"] != lock["lock_sha256"] for row in overlap_rows)
                or any(row["atlas_map_sha256"] != sha256(ROOT / "brain6/results/global/brain6_72_locked.tsv")
                       for row in overlap_rows)
                or any(row["trait_metadata_sha256"] != sha256(ROOT / "results/atlas/traits.tsv")
                       for row in overlap_rows)):
            raise SystemExit("LAVA overlap family does not bind the locked pair selection and source tables")
        matrix_files = []
        for row in overlap_rows:
            matrix_path = ROOT / row["matrix_path"]
            if sha256(matrix_path) != row["matrix_sha256"]:
                raise SystemExit(f"LAVA overlap matrix hash mismatch: {matrix_path}")
            matrix_files.append({"pair_id": row["pair_id"], "path": row["matrix_path"],
                                 "sha256": row["matrix_sha256"]})
        overlap = {"status": "LDSC_INTERCEPT_DERIVED_PAIRWISE_COVARIANCE",
                   "manifest_path": str(overlap_manifest.relative_to(ROOT)),
                   "manifest_sha256": sha256(overlap_manifest),
                   "atlas_map_sha256": overlap_rows[0]["atlas_map_sha256"],
                   "trait_metadata_sha256": overlap_rows[0]["trait_metadata_sha256"],
                   "matrices": matrix_files,
                   "method": "Use each selected pair's stored cross-trait LDSC intercept as the off-diagonal covariance and its two stored univariate LDSC intercepts as diagonal covariances; LAVA converts this matrix to a correlation matrix.",
                   "participant_level_overlap_verified": False,
                   "assume_zero_overlap": False}

    trait_pairs = [(pair["sleep_trait"], pair["disease_trait"]) for pair in lock["pairs"]]
    traits = sorted({trait for pair in trait_pairs for trait in pair})
    n_loci = locus_rows - 1
    univariate_family = len(traits) * n_loci
    pair_locus_family = len(pairs) * n_loci
    document = {
        "schema_version": 1,
        "analysis_id": f"brain6-lava-local-rg-{args.version}",
        "review_status": ("PREDECLARED_WITH_SOURCE_BOUND_OVERLAP; INTERPRETATION_REVIEW_REQUIRED"
                          if overlap else "PREDECLARED_REQUIRES_OVERLAP_INPUT"),
        "scientific_status": "NOT_RUN",
        "selection_is_post_global_screen": True,
        "selection_sha256": selection_actual,
        "pair_lock_sha256": lock["lock_sha256"],
        "pairs": pairs,
        "trait_ids": traits,
        "locus_definition": {
            "path": str(loci.relative_to(ROOT)),
            "sha256": locus_sha,
            "n_loci": n_loci,
            "build": "GRCh37/hg19",
            "rule": "Use all and only the pinned 2,495 loci; do not select loci from Brain6 results.",
        },
        "runtime": {
            "lava_version": "0.1.5",
            "lava_commit": "e729a245f7b6923967a96804fbf5246eadf2d6c6",
            "r_major_minor": "4.3",
            "reference": "official UK Biobank LAVA v1.1",
            "reference_prefix": "ref/lava/ukb_v1.1/lava-ukb-v1.1",
            "reference_must_pass_repository_hash_size_and_index_validation": True,
        },
        "univariate_gate": {
            "unit": "unique_trait_by_locus",
            "n_tests": univariate_family,
            "familywise_alpha": 0.05,
            "p_threshold_strictly_less_than": 0.05 / univariate_family,
            "eligibility_rule": "Both traits in a pair must pass the univariate gate at that locus.",
        },
        "bivariate_family": {
            "unit": "locked_pair_by_locus",
            "n_slots": pair_locus_family,
            "correction": "Benjamini-Hochberg across the fixed 12,475-slot family.",
            "non_tested_slot_p_for_family_correction": 1.0,
            "numerical_failure_p_for_family_correction": 1.0,
            "univariate_gate_miss_p_for_family_correction": 1.0,
            "no_overlap_input_p_for_family_correction": 1.0,
            "report_all_statuses": True,
        },
        "overlap_gate": {
            "status": overlap["status"] if overlap else "REQUIRED_REVIEW_PENDING",
            "pairwise_sample_overlap_file_manifest": overlap,
            "assume_zero_overlap": False,
            "participant_level_overlap_is_verified": False,
            "bivariate_execution_permitted_without_method_specific_overlap_file": False,
        },
        "execution": {
            "one_locus_per_process": True,
            "maximum_locus_failure_fraction": 0.01,
            "maximum_univariate_untested_fraction": 0.05,
            "random_seed": 20260922,
            "reuse_track_b_brain6_local_results": False,
        },
    }
    write_json(args.out, document)
    with atomic_text(checksum_path) as stream:
        stream.write(f"{sha256(args.out)}  {args.out.name}\n")
    print(f"Frozen {len(pairs)} pairs × {n_loci} loci; {len(traits)} unique traits × {n_loci} univariate tests.")
    print(f"family={args.out} sha256={sha256(args.out)}")


if __name__ == "__main__":
    main()
