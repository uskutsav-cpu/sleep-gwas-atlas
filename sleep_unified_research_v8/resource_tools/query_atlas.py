#!/usr/bin/env python3
"""Read the frozen candidate release with exact family/phenotype joins.

No estimation, phenotype substitution, threshold selection, or network access.
The candidate remains subject to source-specific release approval.
"""
import argparse
import csv
import json
from pathlib import Path

TABLES = {"core": "core_rg.tsv", "extension": "extension_rg.tsv",
          "validation": "external_validation_rg.tsv"}


def read_table(root, name):
    with (root / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def unique_index(rows, keys):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        if key in result:
            raise ValueError("Ambiguous record key: " + repr(key))
        result[key] = row
    return result


def query(root, family, sleep_trait, outcome_trait):
    estimates = unique_index(read_table(root, TABLES[family]),
                             ("sleep_trait", "outcome_trait"))
    record = estimates.get((sleep_trait, outcome_trait))
    if record is None:
        return {"status": "NO_ESTIMATE_IN_THIS_RELEASE_TABLE",
                "analysis_family": family, "sleep_trait": sleep_trait,
                "outcome_trait": outcome_trait,
                "interpretation": "This lookup does not establish that the relationship was tested, biologically absent, or absent from other resources."}
    constructs = unique_index(read_table(root, "sleep_constructs.tsv"), ("trait_id",))
    phenotypes = unique_index(read_table(root, "outcome_phenotypes.tsv"),
                              ("analysis_family", "trait_id"))
    sleep = constructs[(sleep_trait,)]
    outcome = phenotypes[(family, outcome_trait)]
    if sleep["source_id"] != record["sleep_source_id"]:
        raise ValueError("Sleep source identity mismatch")
    if outcome["source_id"] != record["outcome_source_id"]:
        raise ValueError("Outcome source identity mismatch")
    evidence = {"status": "ESTIMATE_PRESENT", "record_key": [family, sleep_trait, outcome_trait],
                "estimate": record, "sleep_construct": sleep, "outcome_phenotype": outcome,
                "interval_scope": "Marginal 95%; not a simultaneous interval or a contrast test",
                "causal_interpretation": "Not established by genetic correlation",
                "release_status": "Candidate; source-specific rights and metadata approval pending"}
    hold_path = root / "source_direction_hold.json"
    if not hold_path.exists():
        manifest_path = root / "release_manifest.json"
        if manifest_path.exists() and json.loads(manifest_path.read_text()).get("release_version") == "V8_INTERPRETIVE_QUARANTINE_CANDIDATE":
            raise ValueError("Required V8 source-direction hold file is missing")
        hold_path = Path(__file__).resolve().parents[1] / "tables" / "source_direction_hold.json"
    if hold_path.exists():
        holds = json.loads(hold_path.read_text())
        held = {tuple(r["record_key"]) for r in holds["affected_core_records"]}
        if (family, sleep_trait, outcome_trait) in held:
            evidence["status"] = "HISTORICAL_RECORD_WITH_INTERPRETATION_HOLD"
            evidence["directional_interpretation_status"] = "WITHHOLD_PENDING_SOURCE_EFFECT_DIRECTION_RESOLUTION"
            evidence["directional_interpretation_reason"] = holds["basis"]
            evidence["historical_estimate_scope"] = "Preserved reproduction from historical input; signed disease interpretation is not cleared."
    if family == "validation":
        discovery_id = record["discovery_outcome_trait_id"]
        discovery = unique_index(read_table(root, TABLES["extension"]),
                                 ("sleep_trait", "outcome_trait"))[(sleep_trait, discovery_id)]
        evidence["linked_discovery_estimate"] = discovery
        evidence["sleep_source_reused"] = (discovery["sleep_source_id"] == record["sleep_source_id"])
        evidence["validation_scope"] = "External outcome-side validation; not fully independent two-trait replication"
        evidence["difference_test"] = "Not available: joint sampling covariance is not calibrated"
    evidence["sensitivities"] = [r for r in read_table(root, "sensitivities.tsv")
                                if r["kind"] == "rg" and r["sleep_trait"] == sleep_trait
                                and r["outcome_trait"] == outcome_trait]
    return evidence


def search(root, phrase, family=None):
    words = phrase.casefold().split()
    matches = []
    for row in read_table(root, "outcome_phenotypes.tsv"):
        if family and row["analysis_family"] != family:
            continue
        text = " ".join(row.values()).casefold()
        if all(word in text for word in words):
            matches.append(row)
    return {"status": "METADATA_SEARCH", "query": phrase, "matches": matches,
            "interpretation": "Text matches identify candidate records; they do not establish phenotype equivalence or independent samples. Select the exact family and trait ID for an estimate."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    lookup = sub.add_parser("estimate")
    lookup.add_argument("--family", choices=TABLES, required=True)
    lookup.add_argument("--sleep", required=True)
    lookup.add_argument("--outcome", required=True)
    finder = sub.add_parser("search")
    finder.add_argument("phrase")
    finder.add_argument("--family", choices=TABLES)
    args = parser.parse_args()
    if args.command == "estimate":
        result = query(args.release, args.family, args.sleep, args.outcome)
    else:
        result = search(args.release, args.phrase, args.family)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
