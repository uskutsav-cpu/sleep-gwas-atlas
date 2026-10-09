#!/usr/bin/env python3
"""Read-only ranged-acquisition addendum; no target imports or downloads."""
import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
ROOT=PACKAGE.parent
OUT=HERE/"exact_lipid_ranged_review_v2"
CHUNK=64*1024**2


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        while b:=f.read(65536): h.update(b)
    return h.hexdigest()


def evidence(path):
    path=Path(path)
    return {"path":str(path),"sha256":sha(path),"bytes":path.stat().st_size}


def parse_header(path):
    blocks=Path(path).read_text().replace("\r\n","\n").strip().split("\n\n")
    lines=next(b for b in reversed(blocks) if b.startswith("HTTP/")).splitlines()
    d={"status":lines[0].split()[1]}
    for line in lines[1:]:
        if ":" in line:
            k,v=line.split(":",1);d[k.strip().lower()]=v.strip()
    return d


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--hash-acquired",action="store_true");args=ap.parse_args()
    script=PACKAGE/"scripts/14_acquire_exact_lipids_by_range.py"
    plan_path=PACKAGE/"manifests/exact_lipid_ranged_resource_plan_v2.json"
    plan=json.loads(plan_path.read_text())
    registry=ROOT/"config/public_gwas_sources.tsv"
    with registry.open(newline="") as f:
        sources=[r for r in csv.DictReader(f,delimiter="\t") if r["trait_ids"] in {"hdl","ldl","triglycerides"}]
    head_path=PACKAGE/"logs/exact_lipid_upstream_head_v2.json";head=json.loads(head_path.read_text())
    search_path=PACKAGE/"logs/missing_core_ssd_filename_search_v2.json"
    sequential_path=PACKAGE/"logs/hdl_source_acquisition_receipt_v2.json";seq=json.loads(sequential_path.read_text())
    probe_header=PACKAGE/"logs/lipid_range_probe_headers_v2.txt"
    probe_body=Path("/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/range_probe_v2/hdl_first_1MiB.bin")
    probe=parse_header(probe_header)
    checkpoint=json.loads((ROOT/"discovery_extension/core_checkpoint.json").read_text())
    total=sum(int(s["archive_bytes"]) for s in sources)
    identity=lambda r:(r["trait_ids"],r["download_url"],int(r["archive_bytes"]),r["archive_sha256"])
    checks={
        "registry_exact_frozen_sha":sha(registry)==checkpoint["artifact_hashes_sha256"]["config/public_gwas_sources.tsv"],
        "unique_three_source_traits":len(sources)==3 and {s["trait_ids"] for s in sources}=={"hdl","ldl","triglycerides"},
        "plan_sources_exact_registry":plan["sources"]==sources,
        "HEAD_source_identity_exact":{identity(r) for r in sources}=={(r["trait_id"],r["url"],r["expected_bytes"],r["expected_sha256"]) for r in head["sources"]},
        "plan_script_hash":plan["script_sha256"]==sha(script),
        "plan_registry_hash":plan["source_registry_sha256"]==sha(registry),
        "plan_HEAD_hash":plan["HEAD_receipt_sha256"]==sha(head_path),
        "plan_search_hash":plan["filename_search_receipt_sha256"]==sha(search_path),
        "plan_failed_sequential_receipt_hash":plan["sequential_failure_receipt_sha256"]==sha(sequential_path),
        "plan_probe_header_hash":plan["probe_header_sha256"]==sha(probe_header),
        "plan_probe_body_hash":plan["probe_body_sha256"]==sha(probe_body),
        "probe_206":probe["status"]=="206",
        "probe_exact_range_length":probe["content-range"]=="bytes 0-1048575/2257306867" and probe["content-length"]=="1048576" and probe_body.stat().st_size==1048576,
        "probe_original_HEAD_ETag":probe["etag"]==next(r["headers"]["etag"] for r in head["sources"] if r["trait_id"]=="hdl"),
        "failed_sequential_partial_preserved":Path(seq["path"]).stat().st_size==seq["actual_bytes"]==225443840 and seq["actual_sha256"] is None,
        "plan_network_total":plan["expected_new_network_bytes"]==total,
        "plan_footprint":plan["retention_footprint_bytes_including_chunks_and_previous_partial"]==2*total+seq["actual_bytes"]+probe_body.stat().st_size,
        "plan_SSD_reserve":plan["SSD_free_bytes"]>=plan["retention_footprint_bytes_including_chunks_and_previous_partial"]+5*1024**3,
        "plan_eight_workers_64MiB_chunks":plan["workers"]==8 and plan["chunk_bytes"]==CHUNK and plan["maximum_seconds_per_chunk"]==600,
    }
    intervals={}
    for source in sources:
        n=int(source["archive_bytes"])
        parts=[(i,min(i+CHUNK,n)-1) for i in range(0,n,CHUNK)]
        assert parts[0][0]==0 and parts[-1][1]==n-1
        assert all(b+1==c for (_,b),(c,_) in zip(parts,parts[1:]))
        assert sum(b-a+1 for a,b in parts)==n
        intervals[source["trait_ids"]] = parts
    acquired=[];receipts=[]
    for s in sources:
        p=PACKAGE/"logs"/(s["trait_ids"]+"_ranged_acquisition_receipt_v2.json")
        if not p.exists():
            acquired.append({"trait_id":s["trait_ids"],"receipt_present":False});continue
        receipts.append(p);r=json.loads(p.read_text());n=int(s["archive_bytes"])
        chunk_checks=[]
        etag=next(h["headers"]["etag"] for h in head["sources"] if h["trait_id"]==s["trait_ids"])
        for c in r["chunks"]:
            header=parse_header(Path(c["body_path"]).with_suffix(".headers"))
            chunk_checks.append(c["status"]=="PASS_EXACT_RANGE" and c["curl_returncode"]==0 and
                c["actual_bytes"]==c["expected_bytes"]==c["end"]-c["start"]+1 and
                header["status"]=="206" and header["content-range"]==f"bytes {c['start']}-{c['end']}/{n}" and
                header["content-length"]==str(c["actual_bytes"]) and header["etag"]==etag and
                sha(Path(c["body_path"]).with_suffix(".headers"))==c["header_sha256"])
        file=Path(r["path"])
        item={"trait_id":s["trait_ids"],"receipt_present":True,"receipt_path":str(p),
              "metadata_pass":r["status"]=="EXACT_HISTORICAL_SOURCE_REACQUIRED" and
                    r["expected_sha256"]==r["actual_sha256"]==s["archive_sha256"] and
                    r["expected_bytes"]==r["actual_bytes"]==n and file.stat().st_size==n and
                    r["resource_plan_sha256"]==sha(plan_path) and
                    [(c["start"],c["end"]) for c in r["chunks"]]==intervals[s["trait_ids"]] and all(chunk_checks),
              "chunk_count":len(r["chunks"]),"all_chunk_header_receipts_pass":all(chunk_checks),
              "independent_final_hash_performed":False}
        if args.hash_acquired and item["metadata_pass"]:
            before=file.stat();actual=sha(file);after=file.stat()
            item.update(independent_final_hash_performed=True,actual_sha256=actual,
                        independent_final_hash_matches=actual==s["archive_sha256"],
                        stat_unchanged_during_hash=(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns))
        acquired.append(item)
    aggregate=PACKAGE/"logs/exact_lipid_ranged_acquisition_receipt_v2.json"
    aggregate_checks={}
    if aggregate.exists():
        receipts.append(aggregate)
        a=json.loads(aggregate.read_text())
        aggregate_checks={
            "exact_unique_three_sources":len(a["sources"])==3 and {r["trait_id"] for r in a["sources"]}=={s["trait_ids"] for s in sources},
            "exact_original_sizes_hashes":{(r["trait_id"],r["url"],r["actual_bytes"],r["actual_sha256"]) for r in a["sources"]}=={identity(s) for s in sources},
            "all_source_statuses_success":all(r["status"]=="EXACT_HISTORICAL_SOURCE_REACQUIRED" for r in a["sources"]),
            "original_plan_total":a["resource_plan_sha256"]==sha(plan_path) and a["total_source_bytes"]==total,
            "scope_qualified":a["all_exact_expected_hashes"] is True and a["native_reproduction_completed"] is False,
        }
    complete=bool(aggregate_checks) and all(aggregate_checks.values()) and all(r.get("metadata_pass") for r in acquired)
    result={"reviewed_utc":datetime.now(timezone.utc).isoformat(),"checks":checks,
            "all_preflight_checks_pass":all(checks.values()),"reviewed_script_sha256":sha(script),
            "evidence":[evidence(p) for p in [script,plan_path,registry,head_path,search_path,sequential_path,probe_header,probe_body]+receipts],
            "independent_interval_counts":{k:len(v) for k,v in intervals.items()},"acquisition":acquired,"aggregate_receipt_checks":aggregate_checks,
            "all_three_receipts_complete":complete,"all_three_final_files_independently_rehashed":complete and all(r.get("independent_final_hash_matches") and r.get("stat_unchanged_during_hash") for r in acquired),
            "review_code_sha256":sha(__file__),"large_file_hash_buffer_bytes":65536 if args.hash_acquired else 0,
            "no_downloads_or_native_analyses_launched":True,
            "material_admission_faults":[],
            "operational_limit":"A failed future does not cancel already queued futures; executor context waits. Retention/worker/time caps still apply and failed source is not assembled/admitted.",
            "scientific_limit":"Source byte recovery only; preprocessing/LDSC/family inference/independent biological replication are not reproduced by this script."}
    assert result["all_preflight_checks_pass"]
    result["verdict"]="PASS_CURRENT_RANGED_SOURCE_RECOVERY_PREFLIGHT; "+("THREE_SOURCE_RECEIPTS_COMPLETE" if complete else "ACQUISITION_PENDING_OR_INCOMPLETE")
    lines=["# Exact lipid ranged-transfer review: separate v2 addendum", "",
        "Verdict: **the frozen ranged-transfer design passes source-recovery preflight; "+("all three source receipts are complete" if complete else "acquisition is pending or incomplete")+".** No material source-admission fault was found.", "",
        f"Reviewed script SHA-256: `{sha(script)}`. Frozen resource plan SHA-256: `{sha(plan_path)}`. This reviewer did not import or run the acquisition script and launched no downloads/native analyses. All {len(checks)} independent small-receipt/preflight checks pass.", "",
        "The three selected traits are unique and match the exact frozen source registry. HEAD identities and raw-header SHA-256 are bound into preflight. The preserved 1 MiB HDL probe has HTTP206, exact Content-Range/Content-Length and the original HEAD ETag. The failed sequential partial remains 225,443,840 bytes with no admitted SHA; old source files, partials and sealed v1 outputs are preserved.", "",
        f"The separate frozen plan permits eight curl workers, 64 MiB chunks and 600 seconds per chunk. Expected new source network bytes are {total:,}; retained source/chunk/final/previous-partial/probe footprint is {plan['retention_footprint_bytes_including_chunks_and_previous_partial']:,} bytes, with a 5 GiB reserve against {plan['SSD_free_bytes']:,} SSD bytes free. Internal source bytes are zero. Hash/assembly buffers are 4 MiB; this independent reviewer optionally rehashes final files with 64 KiB reads. Chunk header/receipt overhead fits within the reserve; this does not authorize native-analysis resource use.", "",
        "Each chunk uses If-Match plus a byte range and a per-range max-filesize cap on curl>=8.4. Admission requires successful curl, exact received bytes, HTTP206, exact Content-Range including full source size, matching Content-Length and original ETag. The fresh destination, exclusive JSON creation and prior chunk/header/body checks prevent ordinary rerun overwrite. Every passing chunk has a body/header hash receipt.", "",
        "Independent interval construction produces 34 HDL, 35 LDL and 35 triglycerides chunks with no gaps or overlaps and exact total coverage. Assembly sorts by start, compares the full interval sequence, checks summed bytes, rehashes each chunk against its receipt, then writes via exclusive creation. Final admission requires exact full-file size and the original historical SHA-256. ETag/protocol checks are preconditions; final SHA-256 establishes content identity. Retained chunks support auditability and explain the two-copy storage budget.", "",
        "A failed future does not cancel already queued transfers: ThreadPoolExecutor's context waits for submitted work. Thus failure can continue bounded chunk acquisition until pending workers finish. Each transfer remains capped by chunk bytes/time and eight workers; failure exits before assembly/admission. This is an operational limitation, not a false source-admission pathway. A later plan could cancel unstarted futures for faster failure termination.", "",
        "The permitted SSD temporary-namespace search gap remains explicitly recorded; filename/one-archive searches do not prove all possible source bytes absent. Success would restore upstream archive bytes in a new directory, not missing aliases or completed raw-to-munged/native inference. No biological finding, full 396/1,200-family reproduction, or independent replication follows from these hashes.", "",
        "| Trait | Final receipt review | Independent full-file hash |","|---|---|---|"]
    for r in acquired:
        lines.append(f"| {r['trait_id']} | {'exact source/chunk metadata passes' if r.get('metadata_pass') else 'pending/incomplete'} | {'original SHA matches' if r.get('independent_final_hash_matches') else 'not performed'} |")
    lines += ["", "The corresponding JSON records paths/hashes, independent coverage counts and acquisition state. Re-run `python3 sleep_unified_research_v2/reviews/exact_lipid_ranged_review_v2.py` after receipts arrive; `--hash-acquired` additionally verifies completed full-file bytes. Code/report/JSON are sealed by the matching .sha256 file."]
    Path(str(OUT)+".json").write_text(json.dumps(result,indent=2)+"\n")
    Path(str(OUT)+".md").write_text("\n".join(lines)+"\n")
    fs=[Path(__file__),Path(str(OUT)+".md"),Path(str(OUT)+".json")]
    Path(str(OUT)+".sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in fs))
    print(json.dumps({"verdict":result["verdict"],"checks":len(checks),"report_sha256":sha(Path(str(OUT)+".md")),"receipt_sha256":sha(Path(str(OUT)+".json"))},indent=2))


if __name__=="__main__":main()
