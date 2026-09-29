#!/usr/bin/env python3
"""Reconcile exact-reference partial candidates with Phase 3/30x genotype LD."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
FROZEN = ROOT / "brain6/config/shared_locus_rule_v1.json"
V1 = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v1"
RESCUE = ROOT / "brain6/results/loci/partial_chr17_phase3_rescue_v1"
OUT = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v2"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    rule = json.loads(FROZEN.read_text())
    assert rule["clumping"]["only_exact_snp_chr_position_reference_matches_are_clumped"] is True
    assert rule["missing_reference_variant_policy"].startswith("Retain in candidate-variant audit")
    candidates = {(x["pair_id"], x["SNP"]): x for x in tsv(CANDIDATES)}
    assignments = tsv(V1 / "variant_assignments.tsv")
    rescue = {(x["pair_id"], x["SNP"]): x for x in tsv(RESCUE / "rescue_assignments.tsv")}
    leads = tsv(V1 / "diagnostic_leads.tsv")
    regions = tsv(V1 / "diagnostic_regions.tsv")
    assert len(candidates) == len(assignments) == 2297
    assert len(leads) == 21 and len(regions) == 19 and len(rescue) == 31
    prior = json.loads((V1 / "provenance.json").read_text())
    rescue_receipt = json.loads((RESCUE / "provenance.json").read_text())
    for name, digest in prior["outputs_sha256"].items():
        assert sha(V1 / name) == digest
    for name, digest in rescue_receipt["outputs_sha256"].items():
        assert sha(RESCUE / name) == digest
    lead_keys = {(x["pair_id"], x["lead_SNP"]) for x in leads}
    changes = []
    for row in assignments:
        key = (row["pair_id"], row["SNP"])
        source = candidates[key]
        if row["status"] != "NOT_CLUMPED_NO_EXACT_GENOTYPE":
            assert source["reference_status"] == "EXACT_MATCH"
            continue
        if source["reference_status"] == "NO_REFERENCE_VARIANT":
            assert row["lead_SNP"] == "" and row["r2_to_lead"] == ""
            row["status"] = "NOT_CLUMPED_NO_EXACT_LOCKED_REFERENCE"
            changes.append((key, row["status"]))
            continue
        assert source["reference_status"] == "EXACT_MATCH"
        fixed = rescue[key]
        assert fixed["status"] == "CLUMPED_TO_EXISTING_LEAD"
        assert fixed["phase3_exact_genotype"] == "YES"
        assert float(fixed["max_r2_to_earlier_30x_lead"]) >= 0.1
        assert (row["pair_id"], fixed["best_earlier_lead"]) in lead_keys
        row.update(status="LD_CLUMPED", lead_SNP=fixed["best_earlier_lead"],
                   r2_to_lead=fixed["max_r2_to_earlier_30x_lead"],
                   reference="PHASE3_GRCH37_CANDIDATE_TO_30X_GRCH38_LEAD")
        changes.append((key, row["status"]))
    counts = Counter(x["status"] for x in assignments)
    assert counts == {"LEAD": 21, "LD_CLUMPED": 2252, "NOT_CLUMPED_NO_EXACT_LOCKED_REFERENCE": 24}
    assert len(changes) == 32
    assert all(candidates[(x["pair_id"], x["SNP"])]["reference_status"] == "EXACT_MATCH"
               for x in assignments if x["status"] in {"LEAD", "LD_CLUMPED"})
    assert rescue_receipt["counts"] == {
        "missing_30x_candidates": 31, "phase3_exact_genotypes": 8,
        "unresolved_no_exact_genotype": 23, "clumped_to_existing_lead": 8,
        "new_lead_potential": 0,
    }
    concordance = tsv(RESCUE / "lead_cross_panel_concordance.tsv")
    assert len(concordance) == 3 and all(float(x["cross_panel_r2"]) > 0.99 for x in concordance)
    OUT.mkdir(parents=True, exist_ok=False)
    table = OUT / "variant_assignments.tsv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(assignments[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(assignments)
    for name in ("diagnostic_leads.tsv", "diagnostic_regions.tsv"):
        shutil.copyfile(V1 / name, OUT / name)
    inputs = [CANDIDATES, FROZEN, V1 / "variant_assignments.tsv", V1 / "diagnostic_leads.tsv",
              V1 / "diagnostic_regions.tsv", V1 / "provenance.json",
              RESCUE / "rescue_assignments.tsv", RESCUE / "lead_cross_panel_concordance.tsv",
              RESCUE / "provenance.json", Path(__file__)]
    outputs = [table, OUT / "diagnostic_leads.tsv", OUT / "diagnostic_regions.tsv"]
    receipt = {"analysis_id": "brain6_partial_candidates_1000g_clumping_v2",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "scope": "Partial four-pair candidate diagnostic only; protected Track B absent; no final tiers",
               "interpretation": "Every frozen-rule exact-reference candidate has an independently measured 1000G genotype; 24 without exact locked UKB reference remain non-independent by policy. Eight Phase 3 fallback genotypes clump to existing 30x leads. Lead identities and region boundaries are unchanged.",
               "counts": {"candidates": len(assignments), "exact_reference_genotypes": counts["LEAD"] + counts["LD_CLUMPED"],
                          "excluded_no_locked_reference": counts["NOT_CLUMPED_NO_EXACT_LOCKED_REFERENCE"],
                          "rescued_phase3_genotypes": 8, "diagnostic_leads": len(leads),
                          "diagnostic_regions": len(regions)},
               "inputs_sha256": {str(p): sha(p) for p in inputs},
               "outputs_sha256": {p.name: sha(p) for p in outputs}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt["counts"]), flush=True)


if __name__ == "__main__":
    main()
