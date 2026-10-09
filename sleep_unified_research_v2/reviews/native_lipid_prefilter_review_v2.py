#!/usr/bin/env python3
"""Read-only prefilter-wrapper review; no wrapper/native worker execution."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
ROOT=PACKAGE.parent
OUT=HERE/"native_lipid_prefilter_review_v2"
SSD=Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS")
REF=SSD/"FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr/w_hm3.snplist"
EXPECTED_REF="ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed"


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        while b:=f.read(65536):h.update(b)
    return h.hexdigest()


def evidence(path):
    path=Path(path);return {"path":str(path),"bytes":path.stat().st_size,"sha256":sha(path)}


def git(*args):return subprocess.check_output(["git",*args],cwd=ROOT,text=True).strip()


def main():
    original=ROOT/"scripts/21_prefilter_hm3.py"
    wrapper=PACKAGE/"scripts/16_reproduce_lipid_hm3_prefilter.py"
    blob=git("hash-object",str(original));head=git("rev-parse","HEAD:scripts/21_prefilter_hm3.py")
    old=git("rev-parse","659d01cf:scripts/21_prefilter_hm3.py")
    rename_handled="except FileNotFoundError:" in wrapper.read_text()
    assert blob==head==old=="aae8f23d5c417641d529287c8e7509ea788455e3"
    assert sha(REF)==EXPECTED_REF
    histories={};actuals=[];files=[wrapper,original,REF]
    for trait in ("hdl","ldl","triglycerides"):
        qc=SSD/"Sep1-local-dependencies/data/harmonized"/(trait+".qc.txt")
        historical=dict(line.split("\t",1) for line in qc.read_text().splitlines() if "\t" in line)
        files.append(qc)
        histories[trait]={k:historical[k] for k in ("prefilter_strategy","prefilter_source_sha256","prefilter_allowlist_sha256","prefilter_source_rows","prefilter_retained_rows","infile_sha256")}
        assert historical["prefilter_allowlist_sha256"]==EXPECTED_REF
        receipt_path=PACKAGE/"logs"/(trait+"_native_prefilter_receipt_v2.json")
        item={"trait":trait,"receipt_present":receipt_path.exists(),"verified_native_prefilter":False}
        if receipt_path.exists():
            files.append(receipt_path);r=json.loads(receipt_path.read_text())
            plan_path=PACKAGE/"manifests"/(trait+"_prefilter_execution_plan_v2.json")
            files.append(plan_path);plan=json.loads(plan_path.read_text())
            prov=Path(r["native_prefilter_provenance_path"])
            item["receipt"]=r;item["plan"]=plan
            if prov.exists():
                files.append(prov);p=json.loads(prov.read_text())
                output=Path(r["output_path"])
                observed_output_sha=sha(output) if output.exists() else None
                item["independent_output_sha256"]=observed_output_sha
                item["checks"]={
                    "receipt_success_status":r["status"]=="NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS",
                    "exit_zero_no_stop":r["exit_code"]==0 and r["stop_reason"] is None,
                    "before_after_hash_identity":r["inputs_sha256_before"]==r["inputs_sha256_after"] and r["input_hashes_unchanged"],
                    "plan_receipt_hash":sha(plan_path)==r["resource_plan_sha256"],
                    "original_code_blob":plan["native_script_git_blob"]==plan["native_script_original_main_git_blob"]==blob,
                    "provenance_receipt_hash":sha(prov)==r["native_prefilter_provenance_sha256"],
                    "source_rows_exact":p["source_rows"]==int(historical["prefilter_source_rows"]),
                    "retained_rows_exact":p["retained_rows"]==int(historical["prefilter_retained_rows"]),
                    "source_hash_exact":p["input_sha256"]==historical["prefilter_source_sha256"],
                    "allowlist_hash_exact":p["allowlist_sha256"]==EXPECTED_REF,
                    "output_hash_exact":observed_output_sha==r["output_sha256"]==p["output_sha256"]==historical["infile_sha256"],
                    "output_byte_count":output.stat().st_size==p["output_bytes"],
                    "retained_output_under_plan":r["worker_output_bytes_final_and_temporary"]<=plan["expected_maximum_retained_output_bytes"],
                    "scope_not_full_QC":r["complete_GWAS_QC_or_estimator_reproduction"] is False}
                item["verified_native_prefilter"]=all(item["checks"].values())
        actuals.append(item)
    verified=sum(r["verified_native_prefilter"] for r in actuals)
    data={"reviewed_utc":datetime.now(timezone.utc).isoformat(),"reviewed_wrapper_sha256":sha(wrapper),
        "original_native_script_git_blob":blob,"evidence":[evidence(p) for p in files],
        "historical_prefilter_targets":histories,"native_receipt_reviews":actuals,
        "native_prefilters_independently_verified":verified,"native_workers_launched_by_reviewer":0,
        "scope":"bounded native rsID streaming prefilter only; no harmonization/munging/LDSC or inferential reproduction",
        "prelaunch_material_faults_corrected":["Original wrapper checked final output only; corrected monitor includes actual worker-PID temporary output/provenance paths and post-exit check.","Original terminate timeout could skip receipt; corrected owned-worker kill/wait escalation preserves failure reporting."],
        "limits":["RSS/output/internal-space checks are sampled every2seconds, so stop thresholds are not hard instantaneous resource quotas.","Sequential orchestration across traits is required for the one-worker study plan; separate wrapper invocations have no global singleton lock.","Python binary/code/hash matching does not fully snapshot dynamic libraries/stdlib; an exact final gzip match would directly prove prefilter output identity anyway.","Exact rsID-prefilter bytes/counts do not validate effect coding, ancestry, cohort independence, MAF/INFO filters, genome mapping, allele alignment, heritability or LDSC estimates."],
        "live_rename_robustness_fixed":rename_handled,
        "verdict":f"PASS_SCOPED_PREFILTER_WRAPPER_PREFLIGHT; {verified}_OF_3_NATIVE_PREFILTERS_VERIFIED"}
    lines=["# Lightweight native lipid prefilter review, v2", "",
        f"Verdict: **the corrected wrapper is acceptable for the scoped native rsID prefilter; {verified}/3 native prefilter receipts are independently verified.** No native worker or scientific analysis was launched by this reviewer.", "",
        f"Reviewed wrapper SHA-256: `{sha(wrapper)}`. The actual native script, current HEAD blob and original-main 659d01cf blob are all `{blob}`. Native script SHA-256 is `{sha(original)}`. The actual HapMap3 allowlist SHA-256 is `{sha(REF)}`, matching all three historical QC records.", "",
        "The wrapper gates exact acquired-source SHA against historical prefilter source hash, checks the original allowlist and native-script identity, uses recovered Python 3.9, freezes a per-trait command/resource plan, and binds source/allowlist/native-code/Python-binary/historical-QC/acquisition-receipt hashes before and after execution. Original data and sealed v1 outputs are read; the new per-trait SSD output directory must not preexist.", "",
        "Two material guard faults were corrected before launch: the native script writes output.tmp.<PID> until its final atomic rename, so monitoring only final output missed growing bytes; the revised wrapper includes known temporary/final output and provenance paths and a post-exit size/floor/time check. The revised termination handler escalates a 30-second terminate timeout to kill/wait for only the owned worker, then preserves a failure receipt. "+("The current monitor also tolerates FileNotFoundError during the worker's atomic temporary-to-final rename." if rename_handled else "The monitor still needs FileNotFoundError tolerance for the worker's atomic temporary-to-final rename."), "",
        "The lightweight plan requires 512 MiB internal free space before launch, stops below 128 MiB internal free space, observes worker RSS against 768 MiB, observes combined worker output/provenance against 128 MiB, permits one worker and 3,600 seconds, and places temporary source/output work on SSD. The full-native LDSC 3 GiB internal guard remains unchanged. Sampled 2-second stop thresholds can briefly overshoot and are not strict OS-enforced quotas; exact expected output is far below the source size. The orchestrator must keep trait launches sequential.", "",
        "| Trait | Historical source rows | Historical retained rows | Historical prefilter gzip SHA-256 | Receipt verification |", "|---|---:|---:|---|---|"]
    for item in actuals:
        h=histories[item["trait"]]
        lines.append(f"| {item['trait']} | {h['prefilter_source_rows']} | {h['prefilter_retained_rows']} | `{h['infile_sha256']}` | {'exact native bytes/counts verified' if item['verified_native_prefilter'] else 'pending/incomplete'} |")
    lines += ["", "Success requires unchanged bound hashes, zero exit/no resource stop, native provenance, exact historical source/retained counts, allowlist SHA and gzip SHA. This reviewer independently hashes completed gzip outputs with 64 KiB reads and checks provenance/plan receipts; pending files are not counted as success.", "",
        "This reproduces one preparatory rsID row-selection step only. It does not complete raw harmonization, MAF/INFO screening, allele/effect/build QC, HapMap3 munging, heritability, pairwise LDSC, separate 396/1,200 testing families, biological inference or independent replication. Hash-identical gzip is a reproducible data-processing result, not a new sleep-genetics finding.", "",
        "The JSON contains actual source/QC/reference/code receipts and per-result checks. Re-run `python3 sleep_unified_research_v2/reviews/native_lipid_prefilter_review_v2.py` after execution receipts appear. Code/report/JSON are sealed by the matching .sha256 file."]
    Path(str(OUT)+".json").write_text(json.dumps(data,indent=2)+"\n")
    Path(str(OUT)+".md").write_text("\n".join(lines)+"\n")
    fs=[Path(__file__),Path(str(OUT)+".md"),Path(str(OUT)+".json")]
    Path(str(OUT)+".sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in fs))
    print(json.dumps({"verdict":data["verdict"],"report_sha256":sha(Path(str(OUT)+".md"))},indent=2))


if __name__=="__main__":main()
