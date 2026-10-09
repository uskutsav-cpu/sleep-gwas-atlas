#!/usr/bin/env python3
"""Additive independent frozen-extension printed-log audit, standard library only.

No native LDSC invocation, GWAS data read, existing collator import or prior-review edit.
Writes only extension_log_concordance_v1.json and .tsv adjacent to this script.
"""
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[2]
EXT=ROOT/"discovery_extension"
LOGROOT=Path("/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/logs")
OUT=Path(__file__).with_suffix("")
INPUTS={}
ERRORS=[]
COUNTS=collections.Counter()
ROWS=[]
RTOL=5e-14


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(p):
    INPUTS[str(p)]={"sha256":digest(p),"bytes":p.stat().st_size}
    return p


def load(p):return json.loads(record(p).read_text())


def read(p):
    with record(p).open(newline="") as f:
        return [dict(r,_line=i) for i,r in enumerate(csv.DictReader(f,delimiter="\t"),2)]


def check(family,ident,field,observed,expected,path,line):
    a,b=float(observed),float(expected)
    COUNTS[family+"_numeric_fields_checked"]+=1
    if a!=b:COUNTS[family+"_binary64_nonexact_matches"]+=1
    if not (math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=RTOL,abs_tol=0)):
        ERRORS.append(dict(family=family,identity=ident,field=field,table=observed,log=expected,path=str(path),log_line=line))


def bh(ps):
    assert all(math.isfinite(p) and 0<=p<=1 for p in ps)
    order=sorted(range(len(ps)),key=lambda i:(ps[i],i))
    out=[None]*len(ps);running=1.
    for rank in range(len(ps),0,-1):
        j=order[rank-1];running=min(running,len(ps)*ps[j]/rank);out[j]=running
    return out


def oneline(txt,pat):
    m=re.search(pat,txt,re.M);assert m,pat
    return m


def name(path):return Path(path).name.removesuffix(".sumstats.gz")


def main():
    # Read the bounded recovery receipt first; the much larger tar archive is not opened.
    receipt=load(ROOT/"sleep_unified_research_v1/logs/archived_extension_recovery_receipt_v1.json")
    ledger=read(ROOT/"sleep_unified_research_v1/tables/archived_extension_recovery.tsv")
    bypath={r["path"]:r for r in ledger}
    used=[]
    def log_read(path):
        r=bypath[str(path)]
        record(path)
        actual=INPUTS[str(path)]["sha256"]
        assert actual==r["actual_sha256"]
        assert path.stat().st_size==int(r["bytes"])
        if r["expected_sha256"]:assert actual==r["expected_sha256"]
        used.append(dict(path=str(path),ledger_line=r["_line"],sha256=actual,expected_sha256=r["expected_sha256"],status=r["status"],archive=r["archive_path"],archive_member=r["archive_member"]))
        txt=path.read_text()
        assert "SYNTHETIC SMOKE-TEST OUTPUT" not in txt
        return txt
    hpath=EXT/"results/ldsc/extension_trait_readiness.tsv"
    rgpath=EXT/"results/ldsc/extension_rg_matrix.tsv"
    h2=read(hpath);rg=read(rgpath)
    assert len(h2)==len({r["extension_trait_id"] for r in h2})==100
    assert len(rg)==len({(r["sleep_trait"],r["extension_trait_id"]) for r in rg})==1200
    hs={}
    for row in h2:
        ident=row["extension_trait_id"]
        path=LOGROOT/"h2"/f"h2_{ident}.log"
        txt=log_read(path)
        source=oneline(txt,r"^--h2 (.+) \\\s*$")[1]
        assert name(source)==ident
        hm=oneline(txt,r"^Total (Observed|Liability) scale h2:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)")
        im=oneline(txt,r"^Intercept:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)")
        h,se=float(hm[2]),float(hm[3]);z=h/se
        ninput=int(oneline(txt,r"^Read summary statistics for (\d+) SNPs\.")[1])
        nreg=int(oneline(txt,r"^After merging with regression SNP LD, (\d+) SNPs remain\.")[1])
        expected={"h2":h,"h2_se":se,"h2_z":z,"LDSC_intercept":float(im[1]),"LDSC_intercept_se":float(im[2]),
            "input_snp_count":ninput,"ldsc_regression_snp_count":nreg,
            "lambda_gc":float(oneline(txt,r"^Lambda GC:\s*([-\d.eE+]+)")[1]),
            "mean_chi2":float(oneline(txt,r"^Mean Chi\^2:\s*([-\d.eE+]+)")[1])}
        ratio=re.search(r"^Ratio:\s*([-\d.eE+]+)",txt,re.M)
        if ratio:
            expected["attenuation_ratio"]=float(ratio[1])
            assert row["attenuation_ratio_status"]=="NUMERIC"
        else:
            assert "Ratio < 0 (usually indicates GC correction)." in txt
            assert row["attenuation_ratio"]=="LT_ZERO" and row["attenuation_ratio_status"]=="NEGATIVE_NOT_ESTIMATED"
            COUNTS["h2_negative_ratio_classifications_checked"]+=1
        for field,value in expected.items():check("h2",ident,field,row[field],value,path,txt[:hm.start()].count("\n")+1)
        assert row["scale"]==hm[1].lower()
        gate=z>=4 and float(im[1])<=1.2
        assert (row["primary_rg_eligibility"]=="PRIMARY_PASS")==gate
        assert (row["pass_h2_z_ge_4"]=="True")==(z>=4)
        assert (row["pass_intercept_le_1.2"]=="True")== (float(im[1])<=1.2)
        hs[ident]=dict(input_count=ninput,regression_count=nreg,h2=h,se=se,z=z,intercept=float(im[1]),gate=gate)
        ROWS.append(dict(family="h2",identity=ident,table_line=row["_line"],log_path=str(path),log_line=txt[:hm.start()].count("\n")+1,
            numeric_checks=len(expected),table_input_snp_count=row["input_snp_count"],log_input_snp_count=ninput,
            table_merged_snp_count=row["ldsc_regression_snp_count"],log_merged_snp_count=nreg,
            table_valid_snp_count="NA",log_valid_snp_count="NA",fixed_width_p="NA",scalar_p="NA",table_p="NA",status="PRINTED_LOG_CONCORDANT"))
    observed={}
    width_p_zero=0;width_p_different=0
    sleep_first_h2=[]
    scalar_fields={"rg":"rg","se":"se","z":"z","h2_obs":"extension_h2_observed","h2_obs_se":"extension_h2_observed_se",
                   "h2_int":"extension_h2_intercept","h2_int_se":"extension_h2_intercept_se",
                   "gcov_int":"cross_trait_LDSC_intercept","gcov_int_se":"cross_trait_LDSC_intercept_se"}
    tableby={(r["sleep_trait"],r["extension_trait_id"]):r for r in rg}
    for sleep in sorted({r["sleep_trait"] for r in rg}):
        path=LOGROOT/"rg"/f"rg_{sleep}.log"
        txt=log_read(path);lines=txt.splitlines()
        at=next(i for i,l in enumerate(lines) if l=="Summary of Genetic Correlation Results")
        header=lines[at+1].split();data=[]
        for i in range(at+2,len(lines)):
            if not lines[i].strip():break
            m=re.match(r"^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$",lines[i]);assert m
            fields=[m[1]+".sumstats.gz",m[2]+".sumstats.gz",*m[3].split()]
            assert len(fields)==len(header)
            data.append(dict(zip(header,fields),_line=i+1))
        assert len(data)==100
        chunks=re.split(r"^Computing rg for phenotype \d+/\d+\s*$","\n".join(lines[:at]),flags=re.M)[1:]
        assert len(chunks)==100
        for logrow,chunk in zip(data,chunks):
            pair=(name(logrow["p1"]),name(logrow["p2"]))
            assert pair[0]==sleep and pair[1] in hs and pair not in observed
            row=tableby[pair]
            scalarp=float(oneline(chunk,r"^P:\s*([-\d.eE+]+)")[1])
            ninput=int(oneline(chunk,r"^Read summary statistics for (\d+) SNPs\.")[1])
            nmerge=int(oneline(chunk,r"^After merging with summary statistics, (\d+) SNPs remain\.")[1])
            nvalid=int(oneline(chunk,r"^(\d+) SNPs with valid alleles\.")[1])
            source=oneline(chunk,r"^Reading summary statistics from (.+) \.\.\.")[1]
            assert name(source)==pair[1]
            for field,col in scalar_fields.items():check("rg",pair,col,row[col],logrow[field],path,logrow["_line"])
            check("rg",pair,"p",row["p"],scalarp,path,logrow["_line"])
            for field,value in [("extension_input_snp_count",ninput),("snp_overlap_after_merge",nmerge),("snp_overlap_valid_alleles",nvalid)]:
                check("rg_snp_counts",pair,field,row[field],value,path,logrow["_line"])
            # An extra printed-P representation is preserved, not substituted for scalar P.
            wp=float(logrow["p"])
            assert math.isfinite(wp) and 0<=wp<=1 and abs(wp-scalarp)<=5.000001e-5
            width_p_zero+=wp==0
            width_p_different+=wp!=scalarp
            COUNTS["rg_fixed_width_p_representations_checked"]+=1
            observed[pair]=dict(scalar_p=scalarp,input_count=ninput,merge_count=nmerge,valid_count=nvalid)
            ROWS.append(dict(family="rg",identity="__".join(pair),table_line=row["_line"],log_path=str(path),log_line=logrow["_line"],
                numeric_checks=13,table_input_snp_count=row["extension_input_snp_count"],log_input_snp_count=ninput,
                table_merged_snp_count=row["snp_overlap_after_merge"],log_merged_snp_count=nmerge,
                table_valid_snp_count=row["snp_overlap_valid_alleles"],log_valid_snp_count=nvalid,
                fixed_width_p=logrow["p"],scalar_p=scalarp,table_p=row["p"],status="PRINTED_LOG_CONCORDANT"))
        first=chunks[0]
        hh=oneline(first,r"Heritability of phenotype 1\n[-]+\nTotal Observed scale h2:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)")
        si=float(oneline(first,r"^Intercept:\s*([-\d.eE+]+)")[1])
        sleep_first_h2.append(dict(sleep=sleep,h2=float(hh[1]),se=float(hh[2]),z=float(hh[1])/float(hh[2]),
            intercept=si,gate_pass=float(hh[1])/float(hh[2])>=4 and si<=1.2,
            evidence="FIRST_PAIR_H2_IN_ORIGINAL_EXTENSION_RG_LOG; NOT_STANDALONE_ALL_SNPS_H2"))
    assert set(observed)==set(tableby)
    q=bh([float(r["p"]) for r in rg])
    for row,val in zip(rg,q):check("bh",(row["sleep_trait"],row["extension_trait_id"]),"extension_fdr",row["extension_fdr"],val,rgpath,row["_line"])
    data=dict(scope="Additive independent recovery/hash and printed-statistic concordance; no native estimator run, no GWAS data read, no full-precision historical recovery.",
        status="PASS" if not ERRORS else "FAIL_NUMERICAL_MISMATCH",
        recovery_receipt_snapshot=receipt,relative_roundtrip_tolerance=RTOL,absolute_tolerance=0,
        original_h2_traits_checked=len(h2),original_rg_logs_checked=12,original_rg_pairs_checked=len(observed),
        all_numeric_check_counts=dict(COUNTS),mismatches=ERRORS,
        h2_qc_pass_count=sum(r["gate"] for r in hs.values()),frozen_1200_bh_significant_count=sum(v<.05 for v in q),
        fixed_width_summary_p_zero_count=width_p_zero,fixed_width_summary_p_different_from_scalar_count=width_p_different,
        h2_log_read_count_range=[min(v["input_count"] for v in hs.values()),max(v["input_count"] for v in hs.values())],
        rg_extension_read_count_counts=dict(collections.Counter(v["input_count"] for v in observed.values())),
        cross_stage_read_count_differences=sum(v["input_count"]!=hs[p[1]]["input_count"] for p,v in observed.items()),
        source_count_stage_note="Pinned standalone h2 parses SNP,N,Z with dropna=True; external rg parses SNP,N,Z,A1,A2 with dropna=False before post-merge dropna. Counts are not the same stage; discrepancy alone does not imply source mismatch.",
        source_content_stage_counts_independently_checked_by_this_reviewer=False,
        sleep_first_pair_h2_evidence=sleep_first_h2,
        log_hash_ledger_matches=len(used),log_hash_status_counts=dict(collections.Counter(r["status"] for r in used)),
        log_historical_expected_hash_count=sum(bool(r["expected_sha256"]) for r in used),recovered_log_hashes=used,inputs=INPUTS)
    OUT.with_suffix(".json").write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    with OUT.with_suffix(".tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(ROWS[0]),delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(ROWS)
    print(json.dumps({k:data[k] for k in ["status","original_h2_traits_checked","original_rg_logs_checked","original_rg_pairs_checked","all_numeric_check_counts","mismatches","h2_qc_pass_count","frozen_1200_bh_significant_count","fixed_width_summary_p_zero_count","fixed_width_summary_p_different_from_scalar_count","h2_log_read_count_range","rg_extension_read_count_counts","log_hash_status_counts","log_historical_expected_hash_count"]},indent=2))
    if ERRORS:raise SystemExit(1)


if __name__=="__main__":main()
