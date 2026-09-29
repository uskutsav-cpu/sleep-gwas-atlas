#!/usr/bin/env python3
"""Map all five-track PLACO candidates to frozen LAVA v3 univariate gates."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from audit_candidate_lava_gate_crosswalk_v1 import (
    CANONICAL, FAMILY, LOCI, PROBE, ROOT, digest, read_loci, read_tsv,
)


SOURCE = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"
GROUPS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1"
OUT = ROOT / "brain6/results/loci/five_track_candidate_lava_gate_v1"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite LAVA gate audit: {OUT}")
    receipt = json.loads((SOURCE / "provenance.json").read_text())
    loci_path, variants_path = SOURCE / "candidate_loci.tsv", SOURCE / "candidate_variants.tsv"
    if (receipt["status"] != "PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY"
            or digest(loci_path) != receipt["candidate_loci_sha256"]
            or digest(variants_path) != receipt["candidate_variants_sha256"]):
        raise ValueError("Five-track candidate loci lack a matching source receipt")
    family = json.loads(FAMILY.read_text())
    gate = family["univariate_gate"]
    threshold = float(gate["p_threshold_strictly_less_than"])
    if gate["n_tests"] != 17465 or abs(threshold - 0.05 / 17465) > 1e-18:
        raise ValueError("Frozen LAVA univariate gate changed")
    candidates = read_tsv(loci_path)
    if len(candidates) != receipt["n_candidate_loci"] or len(candidates) != 25:
        raise ValueError("Incomplete five-track candidate family")
    variant_index = {}
    for row in read_tsv(variants_path):
        key = row["locus_id"], row["SNP"]
        if key in variant_index:
            raise ValueError(f"Duplicate five-track candidate variant: {key}")
        variant_index[key] = row
    loci = read_loci()
    canonical_rows = read_tsv(CANONICAL)
    if len(canonical_rows) != 17465:
        raise ValueError("Canonical LAVA family is incomplete")
    canonical = {}
    for row in canonical_rows:
        key = row["phen"], row["locus_id"]
        if key in canonical:
            raise ValueError(f"Duplicate canonical LAVA cell: {key}")
        canonical[key] = row
    probes = {}
    for row in read_tsv(PROBE):
        key = row["pair_id"], row["locus_id"]
        if key in probes:
            raise ValueError(f"Duplicate pre-existing pairwise probe: {key}")
        probes[key] = row
    group_receipt = json.loads((GROUPS / "provenance.json").read_text())
    member_path = GROUPS / "pair_candidate_members.tsv"
    if (group_receipt["status"] != "PASS_GEOGRAPHIC_RECONCILIATION_DIAGNOSTIC_ONLY"
            or digest(member_path) != group_receipt["output_sha256"][member_path.name]):
        raise ValueError("Cross-pair geographic grouping is not receipt bound")
    group_by_candidate = {}
    for row in read_tsv(member_path):
        key = row["candidate_locus_id"]
        if key in group_by_candidate:
            raise ValueError(f"Duplicate geographic candidate membership: {key}")
        group_by_candidate[key] = row["region_group"]
    if set(group_by_candidate) != {row["locus_id"] for row in candidates}:
        raise ValueError("Geographic groups do not cover the exact five-track candidates")

    output = []
    for candidate in candidates:
        locus_id = candidate["locus_id"]
        mapped = set()
        positions = []
        for snp in candidate["lead_variants"].split(";"):
            variant = variant_index.get((locus_id, snp))
            if variant is None or variant["candidate_status"] != "LEAD":
                raise ValueError(f"Candidate lead lacks an exact variant row: {locus_id}:{snp}")
            chrom, pos = int(variant["CHR"]), int(variant["BP"])
            matches = [block for block in loci if int(block["CHR"]) == chrom
                       and int(block["START"]) <= pos <= int(block["STOP"])]
            if len(matches) != 1:
                raise ValueError(f"Lead has {len(matches)} frozen LAVA blocks: {snp}")
            mapped.add(matches[0]["LOC"])
            positions.append(f"{snp}:{chrom}:{pos}")
        if len(mapped) != 1:
            raise ValueError(f"A candidate's independent leads span LAVA blocks: {locus_id}")
        block = next(iter(mapped))
        traits = candidate["pair_id"].split("__")
        if len(traits) != 2:
            raise ValueError(f"Invalid pair ID: {candidate['pair_id']}")
        cells = [canonical[(trait, block)] for trait in traits]
        passes = [cell["status"] == "TESTED" and float(cell["p"]) < threshold
                  for cell in cells]
        probe = probes.get((candidate["pair_id"], block), {})
        output.append({
            "candidate_locus_id": locus_id,
            "geographic_region_group": group_by_candidate[locus_id],
            "pair_id": candidate["pair_id"],
            "lead_variants": candidate["lead_variants"],
            "lead_chr_bp": ";".join(positions),
            "lava_locus_id": block,
            "trait1": traits[0], "trait1_status": cells[0]["status"],
            "trait1_reason": cells[0]["reason"], "trait1_p": cells[0]["p"],
            "trait1_strict_gate_pass": passes[0],
            "trait2": traits[1], "trait2_status": cells[1]["status"],
            "trait2_reason": cells[1]["reason"], "trait2_p": cells[1]["p"],
            "trait2_strict_gate_pass": passes[1],
            "both_strict_gates_pass": all(passes),
            "preexisting_pair_probe_status": probe.get("status", "NO_PROBE"),
            "preexisting_pair_probe_p": probe.get("p", ""),
            "preexisting_pair_probe_q_full_family": probe.get("q_bh_full_12475", ""),
            "interpretation": "DIAGNOSTIC_ONLY_CANONICAL_FAMILY_FAILED_QC",
        })
    OUT.mkdir(parents=True)
    table_path = OUT / "candidate_lava_gate.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    both = [row for row in output if row["both_strict_gates_pass"]]
    result = {
        "analysis_id": "brain6-five-track-candidate-lava-gate-v1",
        "status": "PASS_DIAGNOSTIC_CROSSWALK_NO_PROMOTION",
        "candidate_rows": len(output),
        "both_strict_gates_pass_rows": len(both),
        "both_strict_gates_pass_by_pair": dict(Counter(row["pair_id"] for row in both)),
        "track_b_both_strict_gates_pass_rows": sum(row["both_strict_gates_pass"]
                                                  for row in output if row["pair_id"] == "insomnia__adhd"),
        "unique_pair_locus_strict_gate_pass": sorted({f"{row['pair_id']}:{row['lava_locus_id']}" for row in both}),
        "canonical_lava_status": "FAILED_QC_NOT_PROMOTED",
        "strict_gate_p": threshold,
        "source_sha256": {str(path.relative_to(ROOT)): digest(path) for path in
                          (loci_path, variants_path, SOURCE / "provenance.json", member_path,
                           GROUPS / "provenance.json", LOCI, FAMILY, CANONICAL, PROBE)},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": digest(table_path),
        "interpretation": "Only a descriptive crosswalk. No bivariate LAVA was run, no family QC was repaired, and no candidate was promoted.",
    }
    (OUT / "provenance.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("candidate_rows", "both_strict_gates_pass_rows",
                       "track_b_both_strict_gates_pass_rows", "unique_pair_locus_strict_gate_pass")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
