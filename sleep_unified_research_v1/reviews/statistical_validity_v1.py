#!/usr/bin/env python3
"""Independent, read-only frozen-output arithmetic/log audit; no LDSC rerun.

Writes only statistical_validity_v1.{json,tsv} adjacent to this script.
Uses Python's standard library; imports no repository analytical code.
"""
import collections
import csv
from array import array
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().with_suffix("")
SSD = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
CORE = ROOT / "sleep_unified_research_v1/sources/recovered"
EXT = ROOT / "discovery_extension"
INPUTS = {}
DETAILS = []


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    INPUTS[str(path)] = {"sha256": digest(path), "size_bytes": path.stat().st_size}
    with path.open(newline="") as handle:
        return [dict(row, _line=i) for i, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2)]


def bh(ps):
    assert all(math.isfinite(p) and 0 <= p <= 1 for p in ps)
    order = sorted(range(len(ps)), key=lambda i: (ps[i], i))
    q = [None] * len(ps)
    running = 1.0
    for rank in range(len(ps), 0, -1):
        index = order[rank - 1]
        running = min(running, ps[index] * len(ps) / rank)
        q[index] = running
    return q


def tail(z):
    return math.erfc(abs(z) / math.sqrt(2))


def sigmask(q):
    return [v < .05 for v in q]


def detail(family, path, row, pair, finding, **values):
    DETAILS.append(dict(family=family, file=str(path), line=row, pair_or_trait=pair,
                        finding=finding, values=json.dumps(values, sort_keys=True)))


def log_rg(path):
    INPUTS[str(path)] = {"sha256": digest(path), "size_bytes": path.stat().st_size}
    lines = path.read_text().splitlines()
    assert not any("SYNTHETIC SMOKE-TEST OUTPUT" in line for line in lines)
    at = next(i for i, line in enumerate(lines) if line == "Summary of Genetic Correlation Results")
    header = lines[at + 1].split()
    scalar_p = [line.split(":", 1)[1].strip() for line in lines[:at] if line.startswith("P:")]
    data = []
    for i in range(at + 2, len(lines)):
        line = lines[i]
        if not line.strip():
            break
        m = re.match(r"^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$", line)
        assert m, (path, i + 1, line)
        fields = [m[1] + ".sumstats.gz", m[2] + ".sumstats.gz", *m[3].split()]
        assert len(fields) == len(header)
        data.append(dict(zip(header, fields), _line=i + 1))
    assert len(scalar_p) == len(data)
    for row, p in zip(data, scalar_p):
        row["scalar_p"] = p
    return data


def unname(path):
    return Path(path).name.removesuffix(".sumstats.gz")


def global_family(path, first, second, qfield):
    data = read(path)
    pairs = [(r[first], r[second]) for r in data]
    ps = [float(r["p"]) for r in data]
    q = bh(ps)
    old = [float(r[qfield]) for r in data]
    zq = bh([tail(float(r["z"])) for r in data])
    ratioq = bh([tail(float(r["rg"]) / float(r["se"])) for r in data])
    sm = sigmask(q)
    return data, dict(rows=len(data), distinct_pairs=len(set(pairs)),
        exact_cartesian=set(pairs) == {(a,b) for a in {p[0] for p in pairs} for b in {p[1] for p in pairs}},
        first_trait_count=len({p[0] for p in pairs}), second_trait_count=len({p[1] for p in pairs}),
        bh_max_abs_error=max(abs(a-b) for a,b in zip(q, old)),
        bh_significant_count=sum(sm), printed_z_bh_significant_count=sum(sigmask(zq)),
        printed_rg_over_se_bh_significant_count=sum(sigmask(ratioq)),
        printed_z_bh_label_changes=sum(a != b for a,b in zip(sm, sigmask(zq))),
        printed_rg_over_se_bh_label_changes=sum(a != b for a,b in zip(sm, sigmask(ratioq))),
        zero_p_count=sum(p == 0 for p in ps),
        nonpositive_se_count=sum(float(r["se"]) <= 0 for r in data))


def covariance_diagnostic():
    """Ratio delta method conditional on archived S,V; not a native rerun."""
    path=SSD/"results/tables/ldsc_covariance_pairs.tsv"
    rows=read(path)
    vp=SSD/"results/tables/ldsc_sampling_covariance_1035x1035.tsv.gz"
    INPUTS[str(vp)]={"sha256":digest(vp),"size_bytes":vp.stat().st_size}
    with gzip.open(vp,"rt",newline="") as handle:
        reader=csv.reader(handle,delimiter="\t")
        header=next(reader)
        ids=header[1:]
        vv=[]
        for k,row in enumerate(reader):
            assert row[0]==ids[k] and len(row)==len(ids)+1
            vv.append(array("d",map(float,row[1:])))
    assert len(vv)==len(ids)==len(rows)==1035
    assert ids==[r["element_id"] for r in rows]
    assert all(float(r["genetic_correlation"])==1 for r in rows if r["trait_1"]==r["trait_2"])
    index={s:i for i,s in enumerate(ids)}
    diagonals={r["trait_1"]:float(r["genetic_covariance"]) for r in rows if r["trait_1"]==r["trait_2"]}
    deriv=[]
    scale_errors=[]
    z_errors=[]
    for row in rows:
        i,j=row["trait_1"],row["trait_2"]
        a,b=diagonals[i],diagonals[j]
        oldse=float(row["genetic_correlation_se"])
        expected=float(row["genetic_covariance_se"])/math.sqrt(a*b)
        scale_errors.append(abs(oldse-expected)/expected)
        z_errors.append(abs(float(row["genetic_covariance_z"])-float(row["genetic_correlation_z"])))
        if i==j:
            continue
        rg=float(row["genetic_correlation"])
        positions=[index[f"{i}__{i}"],index[row["element_id"]],index[f"{j}__{j}"]]
        gradient=[-rg/(2*a),1/math.sqrt(a*b),-rg/(2*b)]
        var=sum(gradient[x]*gradient[y]*vv[positions[x]][positions[y]] for x in range(3) for y in range(3))
        assert var>0
        se=math.sqrt(var)
        deriv.append(dict(pair=row["element_id"],line=row["_line"],rg=rg,fixed_scale_se=oldse,
                          delta_method_se_conditional_archived_V=se,se_ratio_delta_over_fixed=se/oldse))
    for row in sorted(deriv,key=lambda x:abs(math.log(x["se_ratio_delta_over_fixed"])),reverse=True)[:10]:
        detail("core_covariance_export",path,row["line"],row["pair"],"rg_se_missing_denominator_uncertainty",**{k:v for k,v in row.items() if k not in {"pair","line"}})
    ratios=[r["se_ratio_delta_over_fixed"] for r in deriv]
    return dict(scope="Delta-method ratio uncertainty conditional on archived 1082-block S,V; not a new jackknife, not an extension discovery/validation covariance estimator.",
        offdiagonal_count=len(deriv),fixed_scale_se_identity_max_relative_error=max(scale_errors),
        covariance_z_equals_exported_correlation_z_max_abs_error=max(z_errors),
        all_diagonal_rg_one=True,diagonal_exported_nonzero_correlation_se_count=sum(float(r["genetic_correlation_se"])>0 for r in rows if r["trait_1"]==r["trait_2"]),
        delta_method_se_over_fixed_scale_se_min=min(ratios),delta_method_se_over_fixed_scale_se_median=statistics.median(ratios),
        delta_method_se_over_fixed_scale_se_max=max(ratios),
        offdiagonal_se_different_over_20percent_count=sum(abs(r-1)>.2 for r in ratios),
        examples=sorted(deriv,key=lambda x:abs(math.log(x["se_ratio_delta_over_fixed"])),reverse=True)[:10])


def main():
    core_path = CORE / "results/tables/rg_matrix.tsv"
    core, cs = global_family(core_path, "sleep_trait", "disease_trait", "fdr")
    ext_path = EXT / "results/ldsc/extension_rg_matrix.tsv"
    ext, es = global_family(ext_path, "sleep_trait", "extension_trait_id", "extension_fdr")
    panel = read(CORE / "config/analysis_panel.tsv")
    h2_path = CORE / "results/tables/h2_summary.tsv"
    h2 = read(h2_path)
    inclusion = read(CORE / "results/tables/phase1_inclusion.tsv")
    hp = {}
    hlogs = {}
    h2_log_mismatch = []
    nonphysical = []
    for row in h2:
        trait = row["trait"]
        path = SSD / row["input_log"]
        INPUTS[str(path)] = {"sha256": digest(path), "size_bytes": path.stat().st_size}
        txt = path.read_text()
        hm = re.search(r"Total (Observed|Liability) scale h2:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)", txt)
        im = re.search(r"^Intercept:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)", txt, re.M)
        hv, se, intercept = float(hm[2]), float(hm[3]), float(im[1])
        hlogs[trait] = dict(h2=hv,se=se,z=hv/se,intercept=intercept,scale=hm[1])
        expect = {"h2":float(f"{hv:.4g}"), "se":float(f"{se:.4g}"), "intercept":float(f"{intercept:.4g}"), "z":round(hv/se,2)}
        for key, val in expect.items():
            if float(row[key]) != val:
                h2_log_mismatch.append(dict(trait=trait, field=key, table=row[key], log=val))
        hp[trait] = hv/se >= 4 and intercept <= 1.2
        if hm[1] == "Liability" and hv > 1:
            nonphysical.append(dict(trait=trait,h2=hv,se=se,verdict=row["verdict"]))
            detail("core",h2_path,row["_line"],trait,"liability_h2_above_one",**nonphysical[-1])
    inc_mismatch = [r["trait_id"] for r in inclusion if (r["include_phase1"].lower() == "true") != hp[r["trait_id"]]]
    log_rows = {}
    for trait in [r["trait_id"] for r in panel if r["domain"] == "sleep"]:
        for path in [SSD/f"results/logs/rg_{trait}.log", SSD/f"results/logs/rg_qc_failed_sensitivity/rg_{trait}.log"]:
            for row in log_rg(path):
                pair = (unname(row["p1"]), unname(row["p2"]))
                assert pair not in log_rows
                log_rows[pair] = (path, row)
    core_log_mismatches = []
    primary = []
    for row in core:
        pair = (row["sleep_trait"],row["disease_trait"])
        path, log = log_rows[pair]
        for column, field in [("rg","rg"),("se","se"),("z","z"),("p","scalar_p"),("h2_int","h2_int"),("gcov_int","gcov_int")]:
            expected = float(f"{float(log[field]):.4g}")
            if not math.isclose(float(row[column]), expected, rel_tol=1e-14, abs_tol=0):
                core_log_mismatches.append(dict(pair=pair,column=column,table=row[column],log=log[field]))
        is_primary = hp[pair[0]] and hp[pair[1]]
        if is_primary:
            primary.append(row)
        assert is_primary == (row["analysis_tier"] == "PRIMARY_PHASE1")
        if float(row["h2_int"]) > 1.2 and is_primary:
            detail("core",core_path,row["_line"],"__".join(pair),"primary_pairwise_outcome_intercept_above_1.2",
                   pairwise_intercept=float(row["h2_int"]), standalone_intercept=hlogs[pair[1]]["intercept"],
                   frozen_fdr=float(row["fdr"]), log=str(path),log_line=log["_line"])
    qp = bh([float(r["p"]) for r in primary])
    cs.update(primary_count=len(primary), sensitivity_count=len(core)-len(primary),
        primary_fdr_max_abs_error=max(abs(float(r["fdr_primary_phase1"])-q) for r,q in zip(primary,qp)),
        primary_significant_count_using_frozen_396=sum(float(r["fdr"])<.05 for r in primary),
        primary_significant_count_using_372=sum(q<.05 for q in qp),
        sensitivity_significant_count_using_frozen_396=sum(float(r["fdr"])<.05 for r in core if r["analysis_tier"]!="PRIMARY_PHASE1"),
        log_rows=len(log_rows),log_numeric_mismatches=core_log_mismatches,
        h2_rows=len(h2),h2_pass_count=sum(hp.values()),h2_fail_traits=[t for t,p in hp.items() if not p],
        h2_log_mismatches=h2_log_mismatch,inclusion_mismatches=inc_mismatch,
        nonphysical_liability_h2=nonphysical,
        primary_pairwise_outcome_intercept_gt1p2_by_trait=dict(collections.Counter(r["disease_trait"] for r in primary if float(r["h2_int"])>1.2)))
    eh_path = EXT / "results/ldsc/extension_trait_readiness.tsv"
    eh = read(eh_path)
    es.update(h2_trait_count=len(eh),h2_gate_pass_count=sum(float(r["h2"])/float(r["h2_se"])>=4 and float(r["LDSC_intercept"])<=1.2 for r in eh),
        standalone_h2_log_present_count=sum((ROOT/r["input_log"]).is_file() for r in eh),
        rg_log_present_count=sum((ROOT/r["input_log"]).is_file() for r in ext),
        pairwise_outcome_intercept_gt1p2=sum(float(r["extension_h2_intercept"])>1.2 for r in ext),
        pairwise_outcome_h2_z_lt4=sum(float(r["extension_h2_observed"])/float(r["extension_h2_observed_se"])<4 for r in ext),
        maximum_abs_cross_trait_intercept=max(abs(float(r["cross_trait_LDSC_intercept"])) for r in ext),
        minimum_valid_snp_overlap=min(int(r["snp_overlap_valid_alleles"]) for r in ext))
    for row in ext:
        if float(row["extension_h2_intercept"])>1.2 or float(row["extension_h2_observed"])/float(row["extension_h2_observed_se"])<4:
            detail("extension",ext_path,row["_line"],row["sleep_trait"]+"__"+row["extension_trait_id"],"pairwise_outcome_qc_warning",
                   pairwise_h2_z=float(row["extension_h2_observed"])/float(row["extension_h2_observed_se"]),
                   pairwise_intercept=float(row["extension_h2_intercept"]),frozen_fdr=float(row["extension_fdr"]))
    rep_path = EXT / "results/replication/replication_results.tsv"
    rep = read(rep_path)
    manifest = read(EXT / "config/replication_manifest.tsv")
    locks={}
    for name in ["replication_candidate_family.lock.json","replication_manifest.lock.json"]:
        path=EXT/"config"/name
        INPUTS[str(path)]={"sha256":digest(path),"size_bytes":path.stat().st_size}
        locks[name]=json.loads(path.read_text())
        assert locks[name]["pair_count"]==217
        assert locks[name]["pair_ids_in_locked_order"]==[r["pair_id"] for r in rep]
    assert locks["replication_manifest.lock.json"]["manifest_sha256"]==digest(EXT/"config/replication_manifest.tsv")
    assert locks["replication_manifest.lock.json"]["candidate_family_lock_sha256"]==digest(EXT/"config/replication_candidate_family.lock.json")
    rgh = read(EXT / "results/replication/replication_source_h2.tsv")
    source_h2_mismatches=[]
    for row in rgh:
        path=ROOT/row["input_log"]
        INPUTS[str(path)]={"sha256":digest(path),"size_bytes":path.stat().st_size}
        txt=path.read_text()
        hm=re.search(r"Total Observed scale h2:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)",txt)
        im=re.search(r"^Intercept:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)",txt,re.M)
        expected={"h2":float(hm[1]),"h2_se":float(hm[2]),"h2_z":float(hm[1])/float(hm[2]),"LDSC_intercept":float(im[1])}
        for field,value in expected.items():
            if not math.isclose(float(row[field]),value,rel_tol=1e-14):
                source_h2_mismatches.append(dict(source=row["replication_source_id"],field=field,table=row[field],log=value))
    replogs = {}
    for path in sorted((EXT/"logs/replication/rg").glob("rg_replication_*.log")):
        for row in log_rg(path):
            key = (unname(row["p1"]), unname(row["p2"]))
            assert key not in replogs
            replogs[key] = (path,row)
    available = [r for r in rep if r["replication_p"] != "NA"]
    threshold = .05/len(rep)
    heterops = []
    rm = []
    class_mismatches = []
    discovery_map = {(r["sleep_trait"],r["extension_trait_id"]):r for r in ext}
    for row in rep:
        pair = (row["sleep_trait"],row["extension_trait_id"])
        d = discovery_map[pair]
        assert float(row["discovery_rg"])==float(d["rg"]) and float(row["discovery_se"])==float(d["se"])
        if row["replication_p"]=="NA":
            continue
        path,log = replogs[(row["sleep_trait"],row["replication_source_id"])]
        for field,col in [("rg","replication_rg"),("se","replication_se"),("z","replication_z"),("gcov_int","cross_trait_LDSC_intercept")]:
            if float(log[field]) != float(row[col]):
                rm.append(dict(pair=row["pair_id"],field=col,table=row[col],log=log[field]))
        p = tail(float(log["z"]))
        if not math.isclose(p,float(row["replication_p"]),rel_tol=1e-13):
            rm.append(dict(pair=row["pair_id"],field="replication_p",table=row["replication_p"],log=p))
        direction = float(row["replication_rg"])*float(row["discovery_rg"]) > 0
        cls = "REPLICATED" if direction and p<threshold else "DIRECTIONALLY_CONCORDANT" if direction else "FAILED_REPLICATION"
        if row["replication_class"] != cls:
            class_mismatches.append(row["pair_id"])
        ds,rs = float(row["discovery_se"]),float(row["replication_se"])
        delta = float(row["replication_rg"])-float(row["discovery_rg"])
        hz = delta/math.sqrt(ds*ds+rs*rs)
        heterop=tail(hz)
        heterops.append(heterop)
        assert math.isclose(heterop,float(row["heterogeneity_p"]),rel_tol=1e-13)
        detail("replication",rep_path,row["_line"],row["pair_id"],"heterogeneity_arithmetic_assumes_covariance_zero",
               delta=delta,se_difference_cov0=math.sqrt(ds*ds+rs*rs),z_difference_cov0=hz,p_difference_cov0=heterop,
               frozen_class=row["replication_class"],positive=row["replication_class"]=="REPLICATED",log=str(path),log_line=log["_line"])
    hfail = [r for r in rgh if r["primary_status"] != "PASS"]
    bysource={r["replication_source_id"]:r for r in rgh}
    qc=collections.Counter("intercept" if float(bysource[r["replication_source_id"]]["LDSC_intercept"])>1.2 else "h2_z" for r in rep if r["replication_class"]=="UNDERPOWERED")
    positives=[r for r in rep if r["replication_class"]=="REPLICATED"]
    rs = dict(rows=len(rep),unique_pair_ids=len({r["pair_id"] for r in rep}),
        both_replication_lock_memberships_and_hashes_match=True,
        manifest_order_match=[r["pair_id"] for r in rep]==[r["pair_id"] for r in manifest],
        classes=dict(collections.Counter(r["replication_class"] for r in rep)),alpha_exact=threshold,
        estimated_pair_count=len(available), log_pair_count=len(replogs),log_numeric_mismatches=rm,class_mismatches=class_mismatches,
        source_h2_count=len(rgh),source_pass_count=len(rgh)-len(hfail),source_fail_count=len(hfail),
        source_h2_log_numeric_mismatches=source_h2_mismatches,
        qc_excluded_pair_reason_counts=dict(qc),
        positive_unique_outcome_count=len({r["extension_trait_id"] for r in positives}),positive_unique_sleep_count=len({r["sleep_trait"] for r in positives}),
        positive_phenotype_match_counts=dict(collections.Counter(r["phenotype_match_status"] for r in manifest if r["pair_id"] in {p["pair_id"] for p in positives})),
        heterogeneity_nominal_count_all_41=sum(p<.05 for p in heterops),heterogeneity_bh_41_count=sum(q<.05 for q in bh(heterops)),
        heterogeneity_bonferroni_41_count=sum(p<.05/len(available) for p in heterops),
        heterogeneity_nominal_count_positive_23=sum(float(r["heterogeneity_p"])<.05 for r in positives))
    checkpoints=json.loads((EXT/"core_checkpoint.json").read_text())
    recovery_hash={p:digest(CORE/p)==expected for p,expected in checkpoints["artifact_hashes_sha256"].items() if (CORE/p).is_file()}
    receipt=dict(scope="Independent standard-library frozen-table and printed-log audit only; no native LDSC rerun or source-level validation.",
        core=cs,extension=es,replication=rs,covariance_export=covariance_diagnostic(),checkpoint_recovered_exact_hash_match=recovery_hash,inputs=INPUTS)
    OUT.with_suffix(".json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    with OUT.with_suffix(".tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=["family","file","line","pair_or_trait","finding","values"],delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(DETAILS)
    print(json.dumps({"core":cs,"extension":es,"replication":rs,"covariance_export":receipt["covariance_export"],"diagnostic_rows":len(DETAILS)},indent=2))


if __name__ == "__main__":
    main()
