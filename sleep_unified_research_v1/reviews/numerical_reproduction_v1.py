#!/usr/bin/env python3
"""Independent native-pilot arithmetic audit. No estimator execution/import.

Writes only numerical_reproduction_v1.json and .tsv beside this script.
"""
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).with_suffix("")
SSD=Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
PILOT=ROOT/"sleep_unified_research_v1/native/core_pilot_v1"
INPUTS={}


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def record(path):
    INPUTS[str(path)]={"sha256":sha(path),"bytes":path.stat().st_size}
    return path


def load(path):
    return json.loads(record(path).read_text())


def tsv(path):
    with record(path).open(newline="") as f:
        return [dict(row,_line=i) for i,row in enumerate(csv.DictReader(f,delimiter="\t"),2)]


def summary(path):
    lines=record(path).read_text().splitlines()
    start=next(i for i,l in enumerate(lines) if l=="Summary of Genetic Correlation Results")
    header=lines[start+1].split()
    ps=[l.split(":",1)[1].strip() for l in lines[:start] if l.startswith("P:")]
    data=[]
    for i in range(start+2,len(lines)):
        if not lines[i].strip():
            break
        m=re.match(r"^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$",lines[i])
        assert m
        row=dict(zip(header,[m[1]+".sumstats.gz",m[2]+".sumstats.gz",*m[3].split()]),_line=i+1)
        data.append(row)
    assert len(ps)==len(data)
    for row,p in zip(data,ps):row["scalar_p"]=p
    return data


def main():
    cap=load(PILOT/"rg_insomnia__bmi.full_precision.json")
    assert len(cap["estimates"])==1
    e=cap["estimates"][0]
    assert e["status"]=="NATIVE_ESTIMATE_RETURNED"
    record(ROOT/"sleep_unified_research_v1/scripts/native_ldsc_capture.py")
    files={}
    deletes={}
    for name in ["hsq1","hsq2","gencov"]:
        paths=list(PILOT.glob("*."+name+".delete"));assert len(paths)==1
        files[name]=str(paths[0])
        deletes[name]=[float(l) for l in record(paths[0]).read_text().splitlines()]
        assert len(deletes[name])==200 and all(math.isfinite(v) for v in deletes[name])
    assert all(v>0 for name in ["hsq1","hsq2"] for v in deletes[name])
    n=len(deletes["gencov"])
    fullratio=e["gencov"]["tot"]/math.sqrt(e["hsq1"]["tot"]*e["hsq2"]["tot"])
    dr=[c/math.sqrt(a*b) for a,b,c in zip(deletes["hsq1"],deletes["hsq2"],deletes["gencov"])]
    pseudo=[n*fullratio-(n-1)*r for r in dr]
    jack=math.fsum(pseudo)/n
    se=math.sqrt(statistics.variance(pseudo)/n)
    z=fullratio/se
    p=math.erfc(abs(z)/math.sqrt(2))
    arithmetic={"rg_ratio":fullratio,"rg_jknife":jack,"rg_se":se,"z":z,"p":p}
    comparisons=[]
    for name,val in arithmetic.items():
        obs=e[name]
        rtol,atol=(1e-12,1e-300) if name=="p" else (1e-12,1e-15)
        good=math.isclose(val,obs,rel_tol=rtol,abs_tol=atol)
        comparisons.append(dict(quantity=name,independent=val,captured=obs,absolute_difference=abs(val-obs),relative_difference=abs(val-obs)/abs(obs),pass_tolerance=good))
        assert good
    hs=[]
    for name in ["hsq1","hsq2","gencov"]:
        tot=e[name]["tot"]
        pv=[n*tot-(n-1)*v for v in deletes[name]]
        computed=math.sqrt(statistics.variance(pv)/n)
        hs.append(dict(quantity=name+".tot_se",independent=computed,captured=e[name]["tot_se"],absolute_difference=abs(computed-e[name]["tot_se"])))
        assert math.isclose(computed,e[name]["tot_se"],rel_tol=1e-12,abs_tol=1e-15)
    core_path=ROOT/"sleep_unified_research_v1/sources/recovered/results/tables/rg_matrix.tsv"
    old=next(r for r in tsv(core_path) if r["sleep_trait"]=="insomnia" and r["disease_trait"]=="bmi")
    old_log=next(r for r in summary(SSD/"results/logs/rg_insomnia.log") if Path(r["p2"]).name=="bmi.sumstats.gz")
    new_log=summary(PILOT/"rg_insomnia__bmi.log")[0]
    log_fields=[k for k in new_log if k not in {"p1","p2","_line"}]
    assert all(float(old_log[k])==float(new_log[k]) for k in log_fields)
    frozen=[]
    mapping={"rg":"rg_ratio","se":"rg_se","z":"z","p":"p"}
    for field,quantity in mapping.items():
        # Original table uses four-significant-digit serialization of already printed logs.
        log_field="scalar_p" if field=="p" else field
        serialized=float(f"{float(old_log[log_field]):.4g}")
        assert float(old[field])==serialized
        frozen.append(dict(quantity=field,frozen=float(old[field]),native=e[quantity],native_minus_frozen=e[quantity]-float(old[field]),historical_printed_log=float(old_log[log_field]),historical_serialization_matches=True))
    comparison=load(ROOT/"sleep_unified_research_v1/manifests/ldsc_native_code_comparison_v1.json")
    code=[]
    for row in comparison:
        a=ROOT.parent/"ldsc-code"/row["file"]
        b=SSD/"ldsc"/row["file"]
        ah,bh=sha(record(a)),sha(record(b))
        assert ah==bh
        code.append(dict(file=row["file"],sha256=ah,independent_byte_identity=True))
    ledger_path=ROOT/"sleep_unified_research_v1/tables/native_input_hash_checks.tsv"
    ledger=tsv(ledger_path)
    input_hashes=[]
    for path in cap["arguments"]["rg"].split(","):
        row=next(r for r in ledger if r["path"]==path and r["kind"]=="core_munged")
        actual=sha(record(Path(path)))
        assert actual==row["actual_sha256"]
        input_hashes.append(dict(path=path,ledger_line=row["_line"],actual_sha256=actual,expected_sha256=row["expected_sha256"],status=row["status"],independent_current_hash_match=True))
    raw=[dict(trait=r["trait_id"],line=r["_line"],kind=r["kind"],status=r["status"],expected_sha256=r["expected_sha256"],actual_sha256=r["actual_sha256"]) for r in ledger if r["trait_id"] in {"insomnia","bmi"} and r["kind"]=="core_raw"]
    assert all(r["expected_sha256"]==r["actual_sha256"] and r["status"]=="MATCH_EXPECTED_SHA256" for r in raw)
    refs=[r for r in ledger if r["kind"]=="official_ld_reference_member" and r["path"].startswith(cap["arguments"]["ref_ld_chr"])]
    refcheck=[]
    # Rehash only the 44 reference members used by this pilot, not dense raw source GWAS.
    for chromosome in range(1,23):
        for suffix in [".l2.ldscore.gz",".l2.M_5_50"]:
            path=Path(cap["arguments"]["ref_ld_chr"])/f"{chromosome}{suffix}"
            row=next(r for r in refs if r["path"]==str(path))
            actual=sha(record(path))
            assert actual==row["actual_sha256"]==row["expected_sha256"]
            refcheck.append(dict(path=str(path),line=row["_line"],sha256=actual,independent_expected_hash_match=True))
    # Preserve source-level incompleteness separately from the admissible processed inputs.
    rawrows=[r for r in ledger if r["kind"]=="core_raw"]
    result=dict(scope="Independent processed-input native-pilot and stock-delete-array audit; no new LDSC execution or source-QC replay.",
        estimator_static_review="capture wrappers call saved stock estimator once, inspect return values, write JSON and return the unchanged result; no estimator-object/argument mutation detected in wrapper",
        instrumented_vs_uninstrumented_native_comparison_performed=False,
        pilot_libraries=cap["libraries"],pilot_python=cap["python"],native_arguments=cap["arguments"],
        independent_arithmetic=comparisons,univariate_and_covariance_delete_se=hs,
        arithmetic_tolerance={"relative":1e-12,"absolute_estimates":1e-15,"absolute_p":1e-300},
        stock_printed_log_identity=True,stock_summary_numeric_field_count=len(log_fields),frozen_table_line=old["_line"],
        frozen_comparison=frozen,original_log_summary_line=old_log["_line"],native_log_summary_line=new_log["_line"],
        code_independent_byte_comparison=code,pilot_processed_input_hash_checks=input_hashes,pilot_raw_hash_ledger_records=raw,
        pilot_raw_sources_rehashed_by_this_reviewer=False,independent_reference_member_hash_checks=refcheck,
        all_raw_ledger_status_counts=dict(collections.Counter(r["status"] for r in rawrows)),
        all_raw_missing_traits=[r["trait_id"] for r in rawrows if r["status"]=="MISSING"],
        block_delete_counts={k:len(v) for k,v in deletes.items()},
        positive_delete_denominators=True,block_boundary_coordinates_captured=False,
        native_ratio_minus_bias_corrected_jackknife=fullratio-jack,
        full_precision_archival_results_available=False,inputs=INPUTS)
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    with OUT.with_suffix(".tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["quantity","independent","captured","absolute_difference","relative_difference","pass_tolerance"],delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(comparisons)
    print(json.dumps({"arithmetic":comparisons,"delete_se":hs,"frozen_comparison":frozen,"code_identical_count":len(code),"reference_members_verified":len(refcheck),"raw_statuses":result["all_raw_ledger_status_counts"],"missing_raw":result["all_raw_missing_traits"]},indent=2))


if __name__=="__main__":main()
