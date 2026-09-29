#!/usr/bin/env python3
"""Crosswalk blocked PLACO candidate regions to frozen LAVA univariate gates.

This is diagnostic only. It never runs bivariate LAVA or changes region tiers.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "brain6/results/loci/candidate_region_decisions_v3_post_b_admission/decisions.tsv"
VARIANTS = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_variants.tsv"
LOCI = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
FAMILY = ROOT / "brain6/config/lava_family_canonical_v3.json"
CANONICAL = (ROOT / "work/lava-canonical-v3-production"
             / "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
             / "results/canonical_family_results.tsv")
PROBE = ROOT / "brain6/results/lava/lava_unique_gate_probe_v1.tsv"
OUTPUT = ROOT / "brain6/results/loci/candidate_lava_gate_crosswalk_v1"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_loci() -> list[dict[str, str]]:
    with LOCI.open(encoding="utf-8") as stream:
        header = stream.readline().split()
        if header != ["LOC", "CHR", "START", "STOP"]:
            raise ValueError("Frozen LAVA locus schema changed")
        loci = [dict(zip(header, line.split(), strict=True)) for line in stream if line.strip()]
    if len(loci) != 2495 or len({r["LOC"] for r in loci}) != len(loci):
        raise ValueError("Frozen LAVA locus family changed")
    return loci


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite crosswalk: {OUTPUT}")
    family = json.loads(FAMILY.read_text())
    gate = family["univariate_gate"]
    threshold = float(gate["p_threshold_strictly_less_than"])
    if gate["n_tests"] != 17465 or abs(threshold - 0.05 / 17465) > 1e-18:
        raise ValueError("Frozen canonical univariate gate changed")

    decisions = read_tsv(DECISIONS)
    if len(decisions) != 19 or any(r["Decision"] != "BLOCKED_EXTERNAL" for r in decisions):
        raise ValueError("Expected the exact 19 blocked post-B candidate regions")
    if len({r["Region"] + "|" + r["Pair"] for r in decisions}) != len(decisions):
        raise ValueError("Duplicate candidate region/pair")
    variants = {}
    for row in read_tsv(VARIANTS):
        key = row["locus_id"], row["SNP"]
        if key in variants:
            raise ValueError(f"Duplicate candidate variant: {key}")
        variants[key] = row
    loci = read_loci()
    canonical_rows = read_tsv(CANONICAL)
    if len(canonical_rows) != 17465:
        raise ValueError("Canonical family is not complete")
    canonical = {}
    for row in canonical_rows:
        key = row["phen"], row["locus_id"]
        if key in canonical:
            raise ValueError(f"Duplicate canonical cell: {key}")
        canonical[key] = row
    probes = {}
    for row in read_tsv(PROBE):
        key = row["pair_id"], row["locus_id"]
        if key in probes:
            raise ValueError(f"Duplicate pre-existing diagnostic pairwise result: {key}")
        probes[key] = row

    output_rows = []
    for region in decisions:
        lead_positions = []
        locus_ids = set()
        for lead in region["lead_variants"].split(";"):
            source = variants.get((region["locus_id"], lead))
            if source is None:
                raise ValueError(f"Missing candidate lead variant: {region['locus_id']}:{lead}")
            chrom, position = int(source["CHR"]), int(source["BP"])
            matching = [locus for locus in loci if int(locus["CHR"]) == chrom
                        and int(locus["START"]) <= position <= int(locus["STOP"])]
            if len(matching) != 1:
                raise ValueError(f"Lead has {len(matching)} matching frozen LAVA blocks: {lead}")
            locus_ids.add(matching[0]["LOC"])
            lead_positions.append(f"{lead}:{chrom}:{position}")
        if len(locus_ids) != 1:
            raise ValueError(f"Candidate region spans multiple lead LAVA blocks: {region['Region']}")
        lava_locus = next(iter(locus_ids))
        traits = region["Pair"].split("__")
        if len(traits) != 2:
            raise ValueError(f"Malformed frozen pair: {region['Pair']}")
        cells = [canonical[(trait, lava_locus)] for trait in traits]
        passes = [cell["status"] == "TESTED" and float(cell["p"]) < threshold
                  for cell in cells]
        probe = probes.get((region["Pair"], lava_locus), {})
        output_rows.append({
            "region": region["Region"], "pair_id": region["Pair"],
            "candidate_locus_id": region["locus_id"],
            "lead_variants": region["lead_variants"],
            "lead_chr_bp": ";".join(lead_positions),
            "lava_locus_id": lava_locus,
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

    OUTPUT.mkdir(parents=True)
    table_path = OUTPUT / "crosswalk.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    unique_gate_pairs = sorted({(row["pair_id"], row["lava_locus_id"])
                                for row in output_rows if row["both_strict_gates_pass"]})
    summary = {
        "analysis_id": "brain6_candidate_lava_gate_crosswalk_v1",
        "scope": "Diagnostic crosswalk only; no bivariate rerun or region promotion",
        "candidate_region_rows": len(output_rows),
        "candidate_lead_rows": sum(len(row["lead_variants"].split(";")) for row in output_rows),
        "candidate_regions_both_strict_gates_pass": sum(row["both_strict_gates_pass"] for row in output_rows),
        "unique_pair_locus_strict_gate_pass": [dict(pair_id=pair, lava_locus_id=loc)
                                               for pair, loc in unique_gate_pairs],
        "pair_counts": dict(Counter(row["pair_id"] for row in output_rows)),
        "canonical_family_status": "FAILED_QC_NOT_PROMOTED",
        "strict_gate_p": threshold,
        "input_sha256": {str(path.relative_to(ROOT)): digest(path)
                         for path in (DECISIONS, VARIANTS, LOCI, FAMILY, CANONICAL, PROBE)},
        "script_sha256": digest(Path(__file__)),
        "crosswalk_sha256": digest(table_path),
    }
    (OUTPUT / "provenance.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: summary[key] for key in
                      ("candidate_region_rows", "candidate_lead_rows",
                       "candidate_regions_both_strict_gates_pass", "unique_pair_locus_strict_gate_pass")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
