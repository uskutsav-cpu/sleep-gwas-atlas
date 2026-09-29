#!/usr/bin/env python3
"""Locate partial PLACO leads in immutable canonical LAVA v3 cells.

This is a diagnostic join only. It cannot promote loci or repair the failed
family-wide LAVA QC gate.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LEADS = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v2/variant_assignments.tsv"
REGIONS = ROOT / "brain6/results/loci/candidate_region_decisions_v2/decisions.tsv"
BLOCKS = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
RUN = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
CELLS = RUN / "results/canonical_family_results.tsv"
DECISION = RUN / "canonical_family_decision.json"
OUT = ROOT / "brain6/results/loci/partial_leads_lava_cells_v1"
FIELDS = ["Region", "Pair", "lead_SNP", "CHR", "BP", "lava_locus_id", "lava_start", "lava_stop",
          "trait1", "trait1_status", "trait1_p", "trait1_reason", "trait2", "trait2_status",
          "trait2_p", "trait2_reason", "both_local_h2_tested", "canonical_family_status"]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def main() -> None:
    decision = json.loads(DECISION.read_text())
    assert decision["overall_status"] == "FAILED_QC_NOT_PROMOTED"
    blocks: dict[int, list[dict[str, str]]] = {}
    with BLOCKS.open(encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter=" ")
        for block in reader:
            blocks.setdefault(int(block["CHR"]), []).append(block)
    cells = {(r["phen"], r["locus_id"]): r for r in tsv(CELLS)}
    region_for_lead = {}
    for region in tsv(REGIONS):
        assert region["Decision"] == "BLOCKED_EXTERNAL"
        for snp in region["lead_variants"].split(";"):
            key = (region["Pair"], snp)
            assert key not in region_for_lead
            region_for_lead[key] = region["Region"]

    rows = []
    for lead in tsv(LEADS):
        if lead["status"] != "LEAD":
            continue
        chrom, bp = int(lead["CHR"]), int(lead["BP"])
        matches = [b for b in blocks[chrom] if int(b["START"]) <= bp <= int(b["STOP"])]
        assert len(matches) == 1, (lead["SNP"], matches)
        block = matches[0]
        trait1, trait2 = lead["pair_id"].split("__")
        c1 = cells[(trait1, block["LOC"])]
        c2 = cells[(trait2, block["LOC"])]
        rows.append({
            "Region": region_for_lead[(lead["pair_id"], lead["SNP"])],
            "Pair": lead["pair_id"], "lead_SNP": lead["SNP"], "CHR": chrom, "BP": bp,
            "lava_locus_id": block["LOC"], "lava_start": block["START"], "lava_stop": block["STOP"],
            "trait1": trait1, "trait1_status": c1["status"], "trait1_p": c1["p"],
            "trait1_reason": c1["reason"], "trait2": trait2, "trait2_status": c2["status"],
            "trait2_p": c2["p"], "trait2_reason": c2["reason"],
            "both_local_h2_tested": c1["status"] == c2["status"] == "TESTED",
            "canonical_family_status": decision["overall_status"],
        })
    assert len(rows) == 21
    assert len({(r["Pair"], r["Region"]) for r in rows}) == 19
    assert set(region_for_lead) == {(r["Pair"], r["lead_SNP"]) for r in rows}
    OUT.mkdir(parents=True, exist_ok=True)
    output = OUT / "lead_lava_cells.tsv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter("BOTH_TESTED" if r["both_local_h2_tested"] else "AT_LEAST_ONE_NOT_RUN" for r in rows)
    pair_region_tested = {}
    for row in rows:
        key = (row["Pair"], row["Region"])
        pair_region_tested[key] = pair_region_tested.get(key, True) and row["both_local_h2_tested"]
    provenance = {
        "analysis_id": "brain6-partial-leads-lava-cells-v1",
        "scope": "Diagnostic mapping of partial-family PLACO leads to immutable LAVA v3 univariate cells",
        "not_a_promotion_gate": True,
        "canonical_family_status": decision["overall_status"],
        "n_leads": len(rows), "n_pair_regions": len({(r["Pair"], r["Region"]) for r in rows}),
        "counts": dict(counts),
        "pair_region_counts": dict(Counter("ALL_LEADS_BOTH_TESTED" if passed else "AT_LEAST_ONE_LEAD_NOT_RUN"
                                           for passed in pair_region_tested.values())),
        "script_sha256": digest(Path(__file__)),
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in (LEADS, REGIONS, BLOCKS, CELLS, DECISION)},
        "output": {"path": str(output.relative_to(ROOT)), "sha256": digest(output)},
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"leads": len(rows), "pair_regions": provenance["n_pair_regions"], "counts": dict(counts)}))


if __name__ == "__main__":
    main()
