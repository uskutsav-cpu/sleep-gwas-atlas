#!/usr/bin/env python3
"""Create the current complete candidate screen without mutating frozen v1 inputs."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
V1 = BASE / "candidate_replacement_screen.tsv"
LITERATURE = BASE / "current_literature_addendum_2026-09-26.tsv"
LITERATURE_PROVENANCE = BASE / "current_literature_addendum_2026-09-26.provenance.json"
OUTPUT = BASE / "candidate_replacement_screen_complete_v2.tsv"
PROVENANCE = BASE / "candidate_replacement_screen_complete_v2.provenance.json"
FIELDS = ("candidate", "PMID", "DOI", "definition", "variant_count", "N", "effective_N",
          "ancestry", "build", "cases_controls", "summary_stats", "fields", "input_processing",
          "license_access", "sample_overlap", "reference_coverage", "decision")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build() -> dict[str, object]:
    if OUTPUT.exists() or PROVENANCE.exists():
        raise FileExistsError("Refusing to overwrite the complete candidate screen v2")
    rows = read_tsv(V1)
    literature = read_tsv(LITERATURE)
    if len(rows) != 5 or len(literature) != 1 or literature[0]["PMID"] != "41922918":
        raise ValueError("Expected five frozen candidates and the dated Portas literature row")
    rows = [{key: row[key] for key in FIELDS} for row in rows]
    dashti = next(row for row in rows if row["candidate"].startswith("Dashti 2019"))
    dashti["decision"] = (
        "NO FULL-SENSITIVITY PROMOTION: completed outcome-blinded 2,495-locus screen adds 6 strict-gate loci, "
        "below the predeclared +250-locus threshold; retain only as a trait-only diagnostic, not replacement or replication"
    )
    source = literature[0]
    rows.append({
        "candidate": source["candidate"], "PMID": source["PMID"], "DOI": source["DOI"],
        "definition": source["definition"], "variant_count": source["variant_count"],
        "N": source["N"], "effective_N": source["effective_N"], "ancestry": source["ancestry"],
        "build": source["build"], "cases_controls": source["cases_controls"],
        "summary_stats": source["summary_stats"], "fields": source["fields"],
        "input_processing": "Not acquired or harmonized; phenotype is not equivalent to self-reported >=9-hour long sleep",
        "license_access": source["access"],
        "sample_overlap": f"Current sleep GWAS: {source['current_sleep_overlap']}; FI GWAS: {source['frailty_overlap']}",
        "reference_coverage": source["reference_coverage"], "decision": source["decision"],
    })
    with OUTPUT.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_complete_candidate_replacement_screen_v2",
        "status": "PASS_COMPLETE_SIX_CANDIDATE_SCREEN",
        "scope": "Outcome-blinded candidate review; preserves the five frozen v1 rows and integrates the dated Portas device-sleep exclusion",
        "decision_rule": "Candidate use is determined by phenotype equivalence, power, access, ancestry/reference compatibility, and documented harmonizability; no downstream sleep-frailty result was used.",
        "inputs_sha256": {
            str(V1.relative_to(ROOT)): sha256(V1),
            str(LITERATURE.relative_to(ROOT)): sha256(LITERATURE),
            str(LITERATURE_PROVENANCE.relative_to(ROOT)): sha256(LITERATURE_PROVENANCE),
        },
        "outputs_sha256": {OUTPUT.name: sha256(OUTPUT)},
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "candidate_count": len(rows),
        "trait_only_candidate_screens_completed": 1,
        "full_sensitivity_candidates_eligible": 0,
        "full_sensitivity_runs_launched": 0,
        "downstream_loci_eligible": 0,
    }
    with PROVENANCE.open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {"status": provenance["status"], "candidate_count": len(rows), "output": str(OUTPUT.relative_to(ROOT))}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
