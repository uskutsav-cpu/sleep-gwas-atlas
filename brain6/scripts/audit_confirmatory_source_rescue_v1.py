#!/usr/bin/env python3
"""Read-only integrity audit for the 2026-09-27 Brain6 source-rescue boundary."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/confirmatory_source_rescue_20260927"
CANON = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
EXPECTED = {
    CANON / "results/canonical_family_results.tsv": "ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352",
    CANON / "canonical_family_decision.json": "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
    ROOT / "brain6/config/lava_family_canonical_v3.json": "f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417",
    ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv": "fffc0e23ab8c111385d78d00168d8dbdbf20b935e7a8fa4b7b6a77e754f3180c",
    BASE / "other_traits/Insomnia_sumstats_Jansenetal.readme_220525.txt": "a5e92e264da88054173fc441e05f4587a6c6134078cc12b69e7294986fcc84a0",
}
SSD_PILOT = Path("/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927")
SSD_EXPECTED = {
    SSD_PILOT / "insomnia-native-n-pilot-v2/materialization.receipt.json": "24ce34356859f1a3f58390c1dc9193c6849cabe2944b1ac947b08d23f53e4f83",
    SSD_PILOT / "insomnia-native-n-pilot-v2/results/exit.json": "4e3c578c526353d3a0ef8108b9e2915f8e2a47c599e0b30f0358b7a7bcfe1aa1",
    SSD_PILOT / "insomnia-native-n-worker2-postmortem-v3/diagnostic.receipt.json": "a023f33abfd63bbdcda1f781b206a61e0f3629672152775970dd099afcbd0457",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    hashes = {}
    for path, expected in EXPECTED.items():
        observed = sha(path)
        if observed != expected:
            raise ValueError(f"Frozen/source hash mismatch: {path}: {observed}")
        hashes[str(path.relative_to(ROOT))] = observed
    ssd_hashes = {}
    for path, expected in SSD_EXPECTED.items():
        observed = sha(path)
        if observed != expected:
            raise ValueError(f"SSD pilot receipt changed: {path}: {observed}")
        ssd_hashes[str(path)] = observed
    exit_receipt = json.loads((SSD_PILOT / "insomnia-native-n-pilot-v2/results/exit.json").read_text())
    diagnostic = json.loads((SSD_PILOT / "insomnia-native-n-worker2-postmortem-v3/diagnostic.receipt.json").read_text())
    if exit_receipt["exit_codes"] != [0, 1, 0, 0] or diagnostic["status"] != "REPRODUCED_FAILURE" \
            or diagnostic["rows"] != 22 or [row["locus_id"] for row in diagnostic["failed"]] != ["950"] \
            or diagnostic["v2_advancement_rule_reopened"] is not False:
        raise ValueError("Insomnia terminal failure/postmortem status changed")

    burden = json.loads((BASE / "phase1/summary.json").read_text())
    traits = rows(BASE / "phase1/trait_failure_burden.tsv")
    if (burden["canonical_tested"], burden["canonical_not_run"], burden["canonical_failed"],
            burden["minimum_required_recovery"], burden["frozen_maximum_not_run"]) != (13745, 3720, 0, 2847, 873):
        raise ValueError("Frozen family burden changed")
    if len(traits) != 7 or sum(int(row["NOT_RUN"]) for row in traits) != 3720:
        raise ValueError("Trait burden does not cover the exact family")
    by_trait = {row["trait_id"]: int(row["NOT_RUN"]) for row in traits}
    if by_trait.get("longsleep") != 1291 or 3720 - sum(
            count for trait, count in by_trait.items() if trait != "longsleep") != 1291:
        raise ValueError("Long-sleep lower bound changed")

    candidates = rows(ROOT / "brain6/results/brain6_bounded_manuscript_v1/candidate_evidence_25.tsv")
    regions = rows(ROOT / "brain6/results/brain6_bounded_manuscript_v1/region_evidence_20.tsv")
    if len(candidates) != 25 or len(regions) != 20:
        raise ValueError("Protected candidate/region row count changed")
    selection = BASE / "other_traits/insomnia_native_n_88_loci.tsv"
    selection_receipt = BASE / "other_traits/insomnia_native_n_88_loci.selection.json"
    pilot_selection = None
    if selection.exists() and selection_receipt.exists():
        chosen = rows(selection)
        receipt = json.loads(selection_receipt.read_text())
        if len(chosen) != 88 or len({row["LOC"] for row in chosen}) != 88 or receipt["table_sha256"] != sha(selection):
            raise ValueError("Prospective insomnia selection receipt mismatch")
        pilot_selection = {"loci": 88, "selection_sha256": sha(selection), "receipt_sha256": sha(selection_receipt)}

    audited_artifacts = [
        BASE / "phase1/summary.json",
        BASE / "phase1/trait_failure_burden.tsv",
        BASE / "phase1/pair_failure_burden.tsv",
        BASE / "phase1/provenance.json",
        BASE / "public_longsleep/SOURCE_EVIDENCE_LEDGER.md",
        BASE / "public_longsleep/source_evidence.tsv",
        BASE / "local_archive/LOCAL_ARCHIVE_AUDIT.md",
        BASE / "local_archive/LOCAL_LONG_SLEEP_CANDIDATES.tsv",
        ROOT / "brain6/results/brain6_bounded_manuscript_v1/candidate_evidence_25.tsv",
        ROOT / "brain6/results/brain6_bounded_manuscript_v1/region_evidence_20.tsv",
    ]
    artifact_hashes = {str(path.relative_to(ROOT)): sha(path) for path in audited_artifacts}

    report = {
        "status": "PASS_READ_ONLY_INTEGRITY",
        "frozen_sha256": hashes,
        "ssd_pilot_receipt_sha256": ssd_hashes,
        "audited_artifact_sha256": artifact_hashes,
        "canonical": {"TESTED": 13745, "NOT_RUN": 3720, "FAILED": 0,
                      "ceiling": 873, "required_recovery": 2847,
                      "unresolved_longsleep_floor_if_others_perfect": 1291},
        "protected_candidates": len(candidates), "reconciled_regions": len(regions),
        "insomnia_selection": pilot_selection,
        "insomnia_v2": {"worker_exit_codes": [0, 1, 0, 0],
                        "postmortem_failed_locus": "950", "advancement": "DENIED"},
        "full_family_gate": "FAILED_QC_NOT_PROMOTED",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
