#!/usr/bin/env python3
"""Audit pinned SuSiE/coloc code and fail closed on missing real locus inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_METHODS = {"SUSIE_RSS_PRIMARY", "COLOC_SUSIE_PRIMARY"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/fine_mapping_sources.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/fine_mapping_colocalization_contract.json"))
    parser.add_argument("--references", type=Path, default=Path("discovery_extension/config/fine_mapping_method_references.tsv"))
    parser.add_argument("--rscript", type=Path, default=Path(".r-env/bin/Rscript"))
    parser.add_argument("--loci", type=Path, default=Path("discovery_extension/results/pleiotropy/novel_shared_loci.tsv"))
    parser.add_argument("--dense-lock", type=Path, default=Path("discovery_extension/provenance/dense_local_inputs.lock.json"))
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_readiness.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/fine_mapping_dependencies.json"))
    args = parser.parse_args()

    sources = read_tsv(args.sources)
    by_method = {row["method_id"]: row for row in sources}
    if set(by_method) != EXPECTED_METHODS or len(sources) != len(EXPECTED_METHODS):
        raise SystemExit("ERROR: fine-mapping source contract must contain exactly the pinned SuSiE-RSS and coloc-SuSiE rows")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    references = read_tsv(args.references)
    required_references = {"SUSIE_MODEL_2020", "SUSIE_RSS_2022", "COLOC_ORIGINAL_2014", "COLOC_PRIORS_2020", "COLOC_SUSIE_2021", "POLYFUN_2020", "FINEMAP_2016"}
    if {row["reference_id"] for row in references} != required_references or len(references) != len(required_references):
        raise SystemExit("ERROR: fine-mapping method reference registry is incomplete or duplicated")
    if contract["fine_mapping"]["primary_method"] != "SUSIE_RSS_PRIMARY" or contract["colocalization"]["primary_method"] != "COLOC_SUSIE_PRIMARY":
        raise SystemExit("ERROR: primary methods differ between source and analysis contracts")
    if contract["fine_mapping"]["software_version"] != by_method["SUSIE_RSS_PRIMARY"]["software_version"]:
        raise SystemExit("ERROR: susieR version differs between source and analysis contracts")
    if contract["colocalization"]["software_version"] != by_method["COLOC_SUSIE_PRIMARY"]["software_version"]:
        raise SystemExit("ERROR: coloc version differs between source and analysis contracts")

    rows: list[dict[str, str]] = []
    all_descriptions_ok = True
    all_archives_ok = True
    for source in sources:
        archive = Path(source["installed_source_archive_path"])
        archive_ok = archive.is_file() and sha256(archive) == source["source_archive_sha256"]
        description = Path(source["installed_description_path"])
        description_ok = description.is_file()
        all_archives_ok = all_archives_ok and archive_ok
        all_descriptions_ok = all_descriptions_ok and description_ok
        rows.append({
            "method_id": source["method_id"],
            "software": source["software"],
            "software_version": source["software_version"],
            "source_commit": source["source_commit"],
            "source_archive_sha256": source["source_archive_sha256"],
            "installed_source_archive_path": str(archive),
            "source_archive_status": "PASS" if archive_ok else "FAIL",
            "installed_description_path": str(description),
            "observed_installed_description_sha256": sha256(description) if description.is_file() else "NA",
            "description_status": "PASS" if description_ok else "FAIL",
            "entrypoints_status": "PENDING",
            "code_status": "PENDING",
            "upstream_loci_status": "PRESENT" if args.loci.is_file() else "BLOCKED_UPSTREAM",
            "dense_input_status": "PENDING",
            "manifest_status": "PRESENT" if args.manifest.is_file() else "NOT_YET_LOCKED",
            "readiness": "PENDING",
        })

    entrypoints_ok = False
    smoke_detail = "Rscript absent"
    if all_archives_ok and all_descriptions_ok and args.rscript.is_file():
        expression = (
            "suppressPackageStartupMessages(library(susieR));"
            "suppressPackageStartupMessages(library(coloc));"
            "ok<-packageVersion('susieR')=='0.14.2' && packageVersion('coloc')=='5.2.3' && "
            "all(vapply(c('susie_rss','susie_get_cs','estimate_s_rss','kriging_rss'),exists,logical(1),where=asNamespace('susieR'),inherits=FALSE)) && "
            "all(vapply(c('coloc.susie','coloc.abf','runsusie','sensitivity'),exists,logical(1),where=asNamespace('coloc'),inherits=FALSE));"
            "cat(ok)"
        )
        smoke = subprocess.run([str(args.rscript), "-e", expression], check=False, text=True, capture_output=True)
        entrypoints_ok = smoke.returncode == 0 and smoke.stdout.strip().endswith("TRUE")
        smoke_detail = (smoke.stdout + smoke.stderr).strip()
    code_ready = all_archives_ok and all_descriptions_ok and entrypoints_ok

    dense_ready = False
    dense_detail = f"dense input lock absent: {args.dense_lock}"
    if args.dense_lock.is_file():
        try:
            dense = json.loads(args.dense_lock.read_text(encoding="utf-8"))
            dense_ready = dense.get("status") == "VERIFIED" and dense.get("full_resolution") is True
            dense_detail = "verified full-resolution dense input lock" if dense_ready else "dense lock is not VERIFIED/full_resolution"
        except (json.JSONDecodeError, OSError) as exc:
            dense_detail = f"invalid dense input lock: {exc}"
    loci_ready = args.loci.is_file()
    manifest_ready = args.manifest.is_file()
    if not code_ready:
        readiness = "BLOCKED_CODE"
    elif not loci_ready or not dense_ready:
        readiness = "CODE_READY_UPSTREAM_LOCUS_AND_DENSE_INPUTS_BLOCKED"
    elif not manifest_ready:
        readiness = "READY_FOR_RESULT_FREE_QTL_CURATION_AND_MANIFEST_LOCK"
    else:
        readiness = "MANIFEST_PRESENT_RUN_LOCK_VALIDATION_REQUIRED"
    for row in rows:
        row["entrypoints_status"] = "PASS" if entrypoints_ok else "FAIL"
        row["code_status"] = "PASS" if code_ready else "FAIL"
        row["dense_input_status"] = "VERIFIED" if dense_ready else "BLOCKED_UPSTREAM"
        row["readiness"] = readiness

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = {
        "schema_version": "1.0.0",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": readiness,
        "analysis_started": False,
        "sources_sha256": sha256(args.sources),
        "contract_sha256": sha256(args.contract),
        "references_sha256": sha256(args.references),
        "runtime_smoke_detail": smoke_detail,
        "code_ready": code_ready,
        "upstream_loci_present": loci_ready,
        "dense_input_verified": dense_ready,
        "dense_input_detail": dense_detail,
        "manifest_present": manifest_ready,
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": "Package readiness is not a fine-mapping or colocalization result; real locus, LD, QTL, and pre-result manifest locks remain mandatory.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FINEMAPPING_PREFLIGHT_{readiness} code_ready={str(code_ready).lower()} dense_ready={str(dense_ready).lower()}")


if __name__ == "__main__":
    main()
