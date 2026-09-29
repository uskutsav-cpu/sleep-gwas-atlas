#!/usr/bin/env python3
"""Summarize sealed prospective insomnia–ADHD PLACO hits as sensitivity only."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.artifacts import verify_current_artifact  # noqa: E402


RUN = Path("/Volumes/Extreme SSD/brain6-work/placo-track-b-prospective-v1")
COLLATED = RUN / "collated"
RESULT = COLLATED / "results.tsv.gz"
STATUS = COLLATED / "status.json"
CONTRACT = RUN / "contract.json"
REGIONS = ROOT / "brain6/results/loci/candidate_region_decisions_v2/decisions.tsv"
OUT = ROOT / "brain6/results/placo/track_b_prospective_sensitivity_v1"
HIT_FIELDS = ["SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "Q_WITHIN_PAIR",
              "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]
REGION_FIELDS = ["Pair", "Region", "candidate_leads", "prospective_b_hit_count",
                 "prospective_b_min_p", "prospective_b_hit_snps", "interpretation"]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_once(path: Path, content: str) -> None:
    if path.exists():
        if path.read_text() != content:
            raise RuntimeError(f"Existing sensitivity artifact differs: {path}")
    else:
        path.write_text(content)


def render_tsv(fields: list[str], rows: list[dict]) -> str:
    import io
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main() -> None:
    receipt = verify_current_artifact(COLLATED)
    contract = json.loads(CONTRACT.read_text())
    status = json.loads(STATUS.read_text())
    if contract["claim_scope"] != "SENSITIVITY_ONLY_NOT_FROZEN_V3_OR_PROTECTED_LEGACY_TRACK_B":
        raise RuntimeError("Prospective result changed scope")
    if status["status"] != "PASS" or receipt["scientific_status"] != "PASS":
        raise RuntimeError("Prospective result did not pass QC")
    if status["headline_threshold"] != 1e-8 or status["total_rows"] != 5_514_402:
        raise RuntimeError("Prospective testing policy or denominator changed")

    regions = []
    with REGIONS.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            chrom, span = row["Region"].split(":")
            start, stop = map(int, span.split("-"))
            if row["Decision"] != "BLOCKED_EXTERNAL":
                raise RuntimeError("Canonical region decision changed")
            regions.append((row, int(chrom.removeprefix("chr")), start, stop))
    if len(regions) != 19:
        raise RuntimeError("Expected 19 partial pair-regions")

    hits = []
    rows_seen = 0
    with gzip.open(RESULT, "rt", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            rows_seen += 1
            if row["headline"] == "True":
                if row["status"] != "TESTED" or float(row["P_PLACO"]) >= 1e-8:
                    raise RuntimeError("Inconsistent prospective headline row")
                hits.append(row)
    if rows_seen != status["total_rows"] or len(hits) != status["headline_variants"]:
        raise RuntimeError("Prospective result and QC denominator differ")

    overlap = []
    for region, chrom, start, stop in regions:
        matched = [hit for hit in hits if int(hit["CHR"]) == chrom and
                   start <= int(hit["BP"]) <= stop]
        overlap.append({
            "Pair": region["Pair"], "Region": region["Region"],
            "candidate_leads": region["lead_variants"],
            "prospective_b_hit_count": len(matched),
            "prospective_b_min_p": min((float(hit["P_PLACO"]) for hit in matched), default="NA"),
            "prospective_b_hit_snps": ";".join(hit["SNP"] for hit in matched),
            "interpretation": "Coordinate overlap only; distinct pair and sensitivity run; no canonical promotion",
        })

    OUT.mkdir(parents=True, exist_ok=True)
    hit_path = OUT / "headline_variants.tsv"
    overlap_path = OUT / "overlap_partial_regions.tsv"
    save_once(hit_path, render_tsv(HIT_FIELDS, hits))
    save_once(overlap_path, render_tsv(REGION_FIELDS, overlap))
    provenance = {
        "analysis_id": "brain6-track-b-prospective-sensitivity-summary-v1",
        "claim_scope": contract["claim_scope"],
        "not_admissible_as_protected_track_b": True,
        "not_a_five_track_family_correction": True,
        "not_a_promotion_gate": True,
        "rows_tested": status["tested"], "headline_variants": len(hits),
        "partial_pair_regions": len(overlap),
        "partial_pair_regions_with_coordinate_overlap": sum(int(row["prospective_b_hit_count"] > 0) for row in overlap),
        "source": {
            "contract_sha256": digest(CONTRACT), "result_sha256": digest(RESULT),
            "receipt_sha256": digest(COLLATED / "receipt.json"),
            "status_sha256": digest(STATUS), "regions_sha256": digest(REGIONS),
        },
        "script_sha256": digest(Path(__file__)),
        "outputs": {path.name: digest(path) for path in (hit_path, overlap_path)},
    }
    save_once(OUT / "provenance.json", json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"headline_variants": len(hits), "pair_regions_with_coordinate_overlap":
                      provenance["partial_pair_regions_with_coordinate_overlap"]}))


if __name__ == "__main__":
    main()
