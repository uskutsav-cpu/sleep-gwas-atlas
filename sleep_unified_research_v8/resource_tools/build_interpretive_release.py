#!/usr/bin/env python3
"""Build a separate interpretive candidate; never write to historical inputs.

Withheld signed fields are absent from the new view. Historical test labels remain
explicitly historical, with no re-estimation or changed multiple-testing family.
"""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--historical-release", type=Path, required=True)
    ap.add_argument("--v8-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    source, root, out = a.historical_release.resolve(), a.v8_root.resolve(), a.out.resolve()
    if source == out or source in out.parents or out in source.parents:
        raise ValueError("Output must be separate from historical release")
    hold = json.loads((root / "tables/source_direction_hold.json").read_text())
    keys = {tuple(r["record_key"]): r for r in hold["affected_core_records"]}
    assert len(keys) == 36
    schema = json.loads((source / "schema.json").read_text())
    before = {name: sha(source / name) for name in schema["files"]}
    out.mkdir(parents=True, exist_ok=True)
    blocked = []
    for name, contract in schema["files"].items():
        if name != "core_rg.tsv":
            shutil.copyfile(source / name, out / name)
            continue
        with (source / name).open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields = reader.fieldnames + ["directional_interpretation_status", "record_use_scope"]
            rows = list(reader)
        assert len(rows) == 396
        for row in rows:
            key = ("core", row["sleep_trait"], row["outcome_trait"])
            row["directional_interpretation_status"] = "QUALIFIED_SOURCE_DEFINED_ESTIMATE"
            row["record_use_scope"] = "Apply original source, QC, marginal uncertainty and testing-family qualifications"
            if key in keys:
                assert row["sleep_source_id"] == keys[key]["sleep_source_id"]
                assert row["outcome_source_id"] == keys[key]["outcome_source_id"]
                for field in ("rg", "z", "ci_lower_95", "ci_upper_95"):
                    row[field] = ""
                row["directional_interpretation_status"] = "WITHHOLD_PENDING_SOURCE_EFFECT_DIRECTION_RESOLUTION"
                row["record_use_scope"] = "Biological correlation interpretation withheld; SE/P/q/threshold/class fields retained only as historical reproduction metadata, not cleared inference"
                row["interpretation"] = row["record_use_scope"]
                blocked.append(key)
        assert set(blocked) == set(keys)
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        contract["fields"] = fields
        contract["held_records"] = 36
        contract["signed_fields_withheld"] = ["rg", "z", "ci_lower_95", "ci_upper_95"]
    result = json.loads((root / "tables/new_exploratory_transport_result.json").read_text())
    assert result["status"] == "NATIVE_COMPLETE_INDEPENDENT_ACTUAL_RESULT_REVIEW_RETAIN_QUALIFIED_EXPLORATORY"
    r = result["result"]
    with (out / "exploratory_source_transport.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(r), lineterminator="\n")
        w.writeheader(); w.writerow(r)
    schema["schema_version"] = "sleep-atlas-interpretive-candidate/2.0"
    schema["files"]["exploratory_source_transport.tsv"] = {"rows": 1, "fields": list(r), "family": "Separate exploratory source/construct transport; no change to historical families"}
    schema["required_interpretation_file"] = "source_direction_hold.json"
    (out / "schema.json").write_text(json.dumps(schema, indent=2) + "\n")
    (out / "source_direction_hold.json").write_text(json.dumps(hold, indent=2) + "\n")
    base_schema = (source / "SCHEMA.md").read_text()
    (out / "SCHEMA.md").write_text(base_schema + "\n## V8 interpretation quarantine\n\nThis is a separate candidate view, not a correction of historical estimates. All 396 core rows remain, but rg, z and both marginal interval bounds are empty for 36 records using the three de Lange bowel sources. The directional status and use-scope fields are mandatory. Their SE, P, q and historical classification fields remain only for provenance; they are not cleared biological inference. A blank signed field is neither zero nor a null result. The original complete 396-pair BH family is unchanged. Historical V4 tables remain immutable.\n\nThe separate one-row exploratory_source_transport.tsv reports FinnGen R13 insomnia with the historically used MVP outcome. Its interval/P are nominal and input-assumption conditional, not independent replication, a calibrated contrast or a new discovery. This file is not part of the historical families.\n\nSource rights, human scientific approval and institutional requirements are still pending; this directory is not a cleared public deposition.\n")
    assert before == {name: sha(source / name) for name in before}
    outputs = {}
    for p in sorted(out.iterdir()):
        if p.is_file() and not p.name.startswith(".") and p.name != "release_manifest.json":
            outputs[p.name] = {"sha256": sha(p), "bytes": p.stat().st_size}
    manifest = {"release_version": "V8_INTERPRETIVE_QUARANTINE_CANDIDATE", "date": "2026-10-11", "status": "NOT_CLEARED_FOR_PUBLIC_DEPOSITION", "historical_input_digests": before, "held_record_count": len(blocked), "held_record_keys": blocked, "historical_inputs_unchanged": True, "historical_families_unchanged": True, "outputs": outputs}
    (out / "release_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": "INTERPRETIVE_VIEW_BUILT", "core_rows": 396, "held": len(blocked), "exploratory_rows": 1, "historical_inputs_unchanged": True}))


if __name__ == "__main__":
    main()
