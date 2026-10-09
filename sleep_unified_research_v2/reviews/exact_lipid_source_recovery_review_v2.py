#!/usr/bin/env python3
"""Independent read-only recovery review; writes only this review's outputs.

Never imports or runs the acquisition script, downloads, or GWAS tools.
--hash-acquired optionally rehashes sealed acquired sources with 64 KiB reads.
"""
import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parent
OUT = HERE / "exact_lipid_source_recovery_review_v2"
TRAITS = {"hdl", "ldl", "triglycerides"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(65536):
            h.update(block)
    return h.hexdigest()


def evidence(path):
    path = Path(path)
    if not path.exists():
        return {"path":str(path), "present":False}
    st = path.stat()
    return {"path":str(path), "present":True, "bytes":st.st_size,
            "sha256":sha(path), "mtime_ns":st.st_mtime_ns}


def identity(row, registry=False):
    return (row["trait_ids"] if registry else row["trait_id"],
            row["download_url"] if registry else row["url"],
            int(row["archive_bytes"] if registry else row["expected_bytes"]),
            row["archive_sha256"] if registry else row["expected_sha256"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hash-acquired", action="store_true")
    args = ap.parse_args()
    registry = ROOT/"config/public_gwas_sources.tsv"
    with registry.open(newline="") as handle:
        selected = [r for r in csv.DictReader(handle, delimiter="\t") if r["trait_ids"] in TRAITS]
    expected = {identity(r, True) for r in selected}
    script = PACKAGE/"scripts/12_recover_exact_lipid_sources.py"
    plan_path = PACKAGE/"manifests/exact_lipid_download_resource_plan_v2.json"
    plan = json.loads(plan_path.read_text())
    head_path = PACKAGE/"logs/exact_lipid_upstream_head_v2.json"
    head = json.loads(head_path.read_text())
    search_path = PACKAGE/"logs/missing_core_ssd_filename_search_v2.json"
    search = json.loads(search_path.read_text())
    checkpoint = json.loads((ROOT/"discovery_extension/core_checkpoint.json").read_text())
    inputs = [script,registry,plan_path,head_path,search_path,
              PACKAGE/"logs/missing_core_archive_search_receipt_v2.json"]
    header_checks = []
    for row in head["sources"]:
        path = PACKAGE/"logs"/(row["trait_id"]+"_upstream_head_v2.txt")
        inputs.append(path)
        header_checks.append({"trait":row["trait_id"], "actual_sha256":sha(path),
                              "expected_sha256":row["header_sha256"],
                              "matches":sha(path)==row["header_sha256"]})
    checks = {
        "current_registry_exact_frozen_hash":sha(registry)==checkpoint["artifact_hashes_sha256"]["config/public_gwas_sources.tsv"],
        "current_registry_three_unique_expected_traits":len(selected)==3 and {r["trait_ids"] for r in selected}==TRAITS and len(expected)==3,
        "all_sources_public_direct_gzip":all(r["access"]=="PUBLIC" and r["archive_member"]=="DIRECT_GZIP" for r in selected),
        "current_HEAD_exact_unique_identity":len(head["sources"])==3 and {identity(r) for r in head["sources"]}==expected,
        "current_HEAD_sizes_returncodes_valid":all(r["returncode"]==0 and r["size_matches"] and r["headers"]["content-length"]==str(r["expected_bytes"]) for r in head["sources"]),
        "raw_HEAD_header_hashes_match":all(r["matches"] for r in header_checks),
        "plan_registry_hash_matches":plan["source_registry_sha256"]==sha(registry),
        "plan_script_hash_matches":plan["script_sha256"]==sha(script),
        "plan_HEAD_hash_matches":plan["HEAD_receipt_sha256"]==sha(head_path),
        "plan_search_hash_matches":plan["filename_search_receipt_sha256"]==sha(search_path),
        "plan_sources_match_registry":plan["sources"]==selected,
        "plan_expected_retained_total_valid":plan["network_bytes_expected"]==plan["maximum_retained_source_bytes"]==sum(int(r["archive_bytes"]) for r in selected),
        "plan_SSD_headroom_valid":plan["SSD_free_bytes_before"]>=plan["maximum_retained_source_bytes"]+plan["SSD_reserve_bytes"],
        "plan_one_worker_no_body_retries":plan["workers"]==1 and plan["body_retry_count"]==0,
        "plan_only_qualified_OS_temp_gap":plan["accepted_search_coverage_gaps"]==search["errors"] and all(r["path"]=="/Volumes/Extreme SSD/.TemporaryItems" and "Operation not permitted" in r["error"] for r in search["errors"]),
        "no_current_lipid_filename_candidates":not any(Path(r["path"]).name in {s["archive_name"] for s in selected}|{s["raw_files"] for s in selected} for r in search["candidates"]),
    }
    # Deterministic predicate counterexample: duplicates still pass set equality
    # when registry and HEAD are both duplicated. Current real data are unique.
    duplicated_selected = [selected[0]]*3
    duplicated_HEAD = [head["sources"][0]]*3
    duplicate_passes_script_identity_predicate = (
        len(duplicated_selected)==3 and len(duplicated_HEAD)==3 and
        {identity(r,True) for r in duplicated_selected}=={identity(r) for r in duplicated_HEAD} and
        all(r["returncode"]==0 and r["size_matches"] and r["headers"]["content-length"]==str(r["expected_bytes"]) for r in duplicated_HEAD))
    acquired = []
    for source in selected:
        receipt_path = PACKAGE/"logs"/(source["trait_ids"]+"_source_acquisition_receipt_v2.json")
        if not receipt_path.exists():
            acquired.append({"trait_id":source["trait_ids"], "receipt_present":False})
            continue
        inputs.append(receipt_path)
        row = json.loads(receipt_path.read_text())
        item = {"trait_id":source["trait_ids"], "receipt_present":True, "receipt":row,
                "exact_receipt_semantics":identity(row)==identity(source,True) and
                row["actual_bytes"]==int(source["archive_bytes"]) and
                row["actual_sha256"]==source["archive_sha256"] and
                row["curl_returncode"]==0 and row["status"]=="EXACT_HISTORICAL_SOURCE_REACQUIRED",
                "independent_large_file_hash_performed":False}
        acquired_path = Path(row["path"])
        item["file_present"] = acquired_path.is_file()
        if item["file_present"]:
            item["current_file_bytes"] = acquired_path.stat().st_size
        header_path = PACKAGE/"logs"/(source["trait_ids"]+"_download_headers_v2.txt")
        if header_path.exists():
            inputs.append(header_path)
            item["download_header_hash_matches"] = sha(header_path)==row["header_sha256"]
        if args.hash_acquired and item["exact_receipt_semantics"] and item["file_present"]:
            # 64 KiB user-space buffer; no source files copied to internal disk.
            before = acquired_path.stat()
            digest = sha(acquired_path)
            after = acquired_path.stat()
            item.update(independent_large_file_hash_performed=True,
                        independent_sha256=digest,
                        independent_hash_matches=digest==source["archive_sha256"],
                        stat_unchanged_during_hash=(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns))
        acquired.append(item)
    aggregate = PACKAGE/"logs/exact_lipid_source_acquisition_receipt_v2.json"
    aggregate_data = json.loads(aggregate.read_text()) if aggregate.exists() else None
    if aggregate.exists():
        inputs.append(aggregate)
    complete = all(r.get("exact_receipt_semantics") for r in acquired) and aggregate_data is not None
    rehashed = complete and all(r.get("independent_hash_matches") and r.get("stat_unchanged_during_hash") for r in acquired)
    result = {"reviewed_utc":datetime.now(timezone.utc).isoformat(),"runtime":sys.version,
              "scope":"read-only source acquisition review; no recovery/native analyses launched",
              "review_code_sha256":sha(__file__),"reviewed_script_sha256":sha(script),
              "checks":checks,"all_small_receipt_checks_pass":all(checks.values()),
              "evidence":[evidence(p) for p in inputs],"HEAD_header_checks":header_checks,
              "acquisition":acquired,"aggregate_receipt":aggregate_data,
              "acquisition_receipts_complete":complete,"all_acquired_files_independently_rehashed":rehashed,
              "large_hash_buffer_bytes":65536 if args.hash_acquired else 0,
              "resource_snapshot":{"internal_free_bytes":shutil.disk_usage("/System/Volumes/Data").free,
                                   "SSD_free_bytes":shutil.disk_usage("/Volumes/Extreme SSD").free},
              "remaining_faults":[{"severity":"P3 hardening; current frozen registry unaffected",
                  "source_lines":"123-125, 137-143","finding":"Three duplicated registry rows and matching duplicated HEAD rows can still pass len/set gate; explicit selected trait set==TRAITS and identity cardinality==3 missing.",
                  "deterministic_predicate_counterexample_passes":duplicate_passes_script_identity_predicate},
                  {"severity":"P3 evidence preservation hardening","source_lines":"175-183",
                   "finding":"Existing source/partial and JSON receipts are preserved, but preexisting download-header text has no existence guard before curl --dump-header. No such preexisting header collision was observed."}],
              "resolved_pre_download_faults":["Empty/wrong HEAD receipt vacuous admission corrected by three-row identity matching and checked header size.","Unenforced retained-size budget corrected by curl>=8.4 requirement and --max-filesize per expected source byte count."],
              "interpretation_limits":["Filename search proves only absence of registered filenames in readable nonexcluded namespaces; no all-SSD/content-identity absence claim.",".TemporaryItems coverage error remains in original receipt and is explicitly qualified in frozen plan.","HEAD length, timestamp and ETag do not prove source byte identity; successful final size and SHA checks do.","Exact source recovery does not reproduce filtering, harmonization, build mapping, munging, LDSC estimates, family correction, cohort independence or new science.","Recovery stores canonical upstream archive names in a new SSD directory; it does not restore missing original raw aliases or overwrite old inputs.","Resource plan bounds downloaded source storage, not all OS swap/log/network overhead or later native-analysis resource needs."],
              "curl_option_reference":"https://curl.se/docs/manpage.html#--max-filesize"}
    assert result["all_small_receipt_checks_pass"]
    status = "exact three-source receipt recovery complete" if complete else "acquisition pending or incomplete; no three-source recovery claim yet"
    result["verdict"] = "SOURCE_RECOVERY_GATES_ACCEPTABLE_FOR_CURRENT_PINNED_INPUTS; "+status.upper()
    body = ["# Exact lipid source recovery: bounded independent technical review, v2", "",
            "Verdict: **the revised gates and frozen resource plan are adequate for the current three pinned source identities; "+status+".**", "",
            f"Reviewed acquisition script SHA-256: `{sha(script)}`. Frozen resource plan SHA-256: `{sha(plan_path)}`. No acquisition/native function was executed by this reviewer. The Python review uses independent 64 KiB-buffer hashing and writes only its own files in `sleep_unified_research_v2/reviews/`.", "",
            f"All {len(checks)} metadata/receipt/registry checks pass. The registry is byte-identical to the frozen checkpoint (`{sha(registry)}`); it contains exactly HDL, LDL and log-triglycerides, all PUBLIC DIRECT_GZIP sources. Current HEAD receipts contain three distinct matching trait/URL/size/SHA identities, successful curl return codes and exact Content-Length values. Their raw header hashes were checked independently. HEAD is a preflight gate; it is not a content-identity test.", "",
            f"The frozen plan permits one worker, no body retries, {plan['maximum_retained_source_bytes']:,} total source bytes, and a {plan['SSD_reserve_bytes']:,}-byte SSD reserve. SSD free space before acquisition was {plan['SSD_free_bytes_before']:,} bytes. The source bodies go directly to a new SSD directory; the script hashes with 4 MiB reads. Actual curl is 8.7.1, and the revised command enforces a per-source --max-filesize cap. curl documents in-transfer size-limit enforcement from 8.4.0; both the installed local manual and [official curl manual](https://curl.se/docs/manpage.html#--max-filesize) support the version gate. This is a source-storage plan, not clearance for later preprocessing/LDSC memory or internal-disk use.", "",
            "Two material faults identified in the first inspected script were corrected before acquisition: an empty or unrelated HEAD list could pass the original all() test, and expected retained bytes were not enforced during GET. The frozen current plan binds the revised script, registry, original HEAD receipt and original filename-search receipt by SHA-256. The earlier HEAD/search receipts retain their older script SHA, which truthfully identifies their generating version.", "",
            "The original search receipt is preserved with its inaccessible `/Volumes/Extreme SSD/.TemporaryItems` entry. The amended gate accepts only that path and the recorded permission error; other errors and any lipid filename candidate still stop acquisition. The new plan explicitly qualifies the inaccessible namespace. There are no lipid filename candidates in the scanned readable namespaces. A separate named-archive receipt reports no target candidates. Neither proves that no differently named source bytes exist elsewhere.", "",
            "Original source locations and frozen v1 outputs are only read. JSON evidence uses exclusive creation, and existing destination source/partial files cause an abort. Successful files are admitted only when curl returns zero, actual bytes exactly equal the expected count, and streaming SHA-256 equals the frozen source hash; a mismatch is preserved and blocks downstream use. The destination retains canonical archive filenames in a new directory rather than recreating old raw aliases.", "",
            "Residual hardening items do not invalidate the current acquisition: (1) the generalized length/set gate still allows three duplicate registry and corresponding HEAD rows if both are changed together; an isolated predicate counterexample confirms this, but the real frozen registry/HEAD are unique; (2) download-header text lacks a preexisting-file guard before --dump-header, unlike source/partial/JSON files. No header collision was observed. Explicit uniqueness assertions and an exclusive header-path check would make those preservation promises hold in these edge cases.", "",
            "| Trait | Acquisition receipt | Independent acquired-file hash |", "|---|---|---|"]
    for item in acquired:
        receipt_status = "exact bytes/SHA receipt verified" if item.get("exact_receipt_semantics") else ("failed receipt preserved" if item["receipt_present"] else "pending")
        hash_status = "SHA matched with unchanged file stat" if item.get("independent_hash_matches") else "not performed"
        body.append(f"| {item['trait_id']} | {receipt_status} | {hash_status} |")
    body += ["", "Acquisition receipt success can support exact historical upstream byte recovery. It does not demonstrate a raw-to-harmonized-to-munged rerun, allele/build/coding correctness, native LDSC estimates, all 396/1,200 results, biological validity or independent replication. No scientific reproduction is inferred from source hashes.", "",
             "Run this read-only audit with `python3 sleep_unified_research_v2/reviews/exact_lipid_source_recovery_review_v2.py`; add `--hash-acquired` only after final source receipts exist and resource headroom permits sequential source reads. The JSON contains file paths/hashes, deterministic gate counterexample and receipt details. The code/report/receipt are sealed by the matching .sha256 file."]
    Path(str(OUT)+".json").write_text(json.dumps(result,indent=2)+"\n")
    Path(str(OUT)+".md").write_text("\n".join(body)+"\n")
    files = [Path(__file__),Path(str(OUT)+".json"),Path(str(OUT)+".md")]
    Path(str(OUT)+".sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"verdict":result["verdict"],"checks":len(checks),
                      "report_sha256":sha(Path(str(OUT)+".md")),
                      "receipt_sha256":sha(Path(str(OUT)+".json"))},indent=2))


if __name__=="__main__":
    main()
