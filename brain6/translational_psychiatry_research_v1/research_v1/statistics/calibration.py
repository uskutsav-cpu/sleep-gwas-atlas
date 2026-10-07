"""Source-specific method diagnostics. Simulations never count as GWAS validation."""
from __future__ import annotations
import ast, csv, datetime, hashlib, json, math, pathlib, types
from collections import OrderedDict
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE_SHA = "98eb8aa9384d880ded4965d3f920844a02f317304408dd884037d6256e0a90df"

def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def write_tsv(name, rows):
    if not rows: return
    with (HERE/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t"); w.writeheader(); w.writerows(rows)

def wilson(k,n):
    z=norm.ppf(.975); p=k/n; den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half

def exact_ld(m,a):
    r=a**np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])
    d,v=np.linalg.eigh(r); order=np.argsort(d)[::-1]
    return r,d[order],v[:,order],(r*r).sum(axis=1)

def upstream_callable():
    p=HERE/"source"/"calculate.py"
    if sha(p)!=SOURCE_SHA: raise ValueError("pinned source hash mismatch")
    tree=ast.parse(p.read_text())
    funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {"nearest_Corr","calLocalCov"}]
    env={"np":np,"pd":pd,"linalg":np.linalg,"sqrt":math.sqrt,"norm":norm,"OrderedDict":OrderedDict,"ld":types.SimpleNamespace(getBlockLefts=lambda coords,dist:np.zeros(len(coords),dtype=int))}
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(p),"exec"),env)
    return env["calLocalCov"]

def reconstructed(t1,t2,d,n1,n2,h1_context,h2_context,meanld,c=0.,cvar=0.):
    """Reconstruct archived calLocalCov literally, with batched eigen-transformed Z."""
    m=len(d); b=len(t1); rootn=math.sqrt(n1*n2)
    h1=(np.mean(t1*t1,axis=1)-1)*m/(n1*meanld)
    h2=(np.mean(t2*t2,axis=1)-1)*m/(n2*meanld)
    rho=(np.sum(t1*t2,axis=1)-c*m)/(meanld*rootn)
    y=t1*t2-c*d
    q=(h1_context[:,None]*d/m+1/n1)*(h2_context[:,None]*d/m+1/n2)
    v1=np.cumsum(d*d/q,axis=1)
    v2=np.cumsum(y/rootn/q,axis=1)
    v3=np.cumsum(y*y/(n1*n2*q*d*d),axis=1)
    base=np.sum(d>1.)
    if base<2: raise ValueError("illustrative LD needs >=2 eigenvalues >1 for source formula")
    counts=np.arange(base,m+1)
    emp=(v3[:,base-1:]-v2[:,base-1:]**2/v1[:,base-1:])/(v1[:,base-1:]*(counts-1))
    theo=1/v1[:,base-1:]
    vmax=np.maximum(emp,theo)
    which=np.argmin(vmax,axis=1)
    variance=m*m*vmax[np.arange(b),which]
    # Archived source uses [:base+min_idx-1], one fewer mode than its variance selection.
    cut=base+which-1
    v4=np.cumsum(d/q,axis=1)[np.arange(b),cut-1]/v1[np.arange(b),cut-1]
    variance+=cvar/(n1*n2)*m*m*v4*v4
    with np.errstate(invalid="ignore",divide="ignore"):
        corr=np.where((h1<0)|(h2<0),np.nan,rho/np.sqrt(h1*h2))
        p=2*norm.sf(np.abs(rho/np.sqrt(variance)))
    return {"rho":rho,"h2_1":h1,"h2_2":h2,"corr":corr,"var":variance,"p":p,"selected_modes_for_variance":base+which,"selected_modes_for_overlap_variance":cut}

def native_fixture(a=.6,h=.001,c=.2,seed=1,n1=386533,n2=225534):
    """Execute unchanged upstream numerical function; genotype/LD API supplies known synthetic LD."""
    rng=np.random.default_rng(seed); m=160; r,d,v,l2=exact_ld(m,a)
    av=d+n1*h*d*d/m; bv=d+n2*h*d*d/m; cross=c*d
    tx=np.sqrt(av)[None,:]*rng.standard_normal((4,m))
    ty=cross/av*tx+np.sqrt(bv-cross*cross/av)[None,:]*rng.standard_normal((4,m))
    zx=(tx@v.T).ravel(); zy=(ty@v.T).ravel()
    g=pd.DataFrame({"CHR":np.ones(4*m),"Z_x":zx,"Z_y":zy})
    ldscore=pd.DataFrame({"L2":np.tile(l2,4)})
    geno=types.SimpleNamespace(ldCorrVarBlocks=lambda left,idx:(l2,r.copy()))
    source=upstream_callable()(0,pd.DataFrame([[1,1,m]]),geno,np.arange(4*m)*.00001,np.arange(1,4*m+1),g,ldscore,n1,n2,c,0.)
    hc1=(np.mean(tx*tx)-1)*m/(n1*np.mean(l2))
    hc2=(np.mean(ty*ty)-1)*m/(n2*np.mean(l2))
    independently=reconstructed(tx[:1],ty[:1],d,n1,n2,np.array([hc1]),np.array([hc2]),np.mean(l2),c)
    errors={k:abs(float(source.iloc[0][k])-float(independently[k][0])) for k in ["rho","h2_1","h2_2","corr","var","p"]}
    return source,independently,errors,(tx,ty,d,np.mean(l2),hc1,hc2)

def freeze():
    paths=[HERE/"calibration_protocol.json",HERE/"calibration.py",HERE/"source"/"calculate.py",HERE/"source"/"pheno.py",HERE/"source"/"heritability.py"]
    if (HERE/"calibration_freeze.json").exists(): raise ValueError("Freeze already exists; do not rewrite")
    record={"frozen_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"no_simulation_outcomes_read":True,"file_sha256":{str(p.relative_to(HERE)):sha(p) for p in paths}}
    (HERE/"calibration_freeze.json").write_text(json.dumps(record,indent=2)+"\n")

def verify_freeze():
    record=json.loads((HERE/"calibration_freeze.json").read_text())
    for p,h in record["file_sha256"].items():
        if sha(HERE/p)!=h: raise ValueError("Frozen code/settings/source changed: "+p)
    return record

def run():
    freeze_record=verify_freeze(); cfg=json.loads((HERE/"calibration_protocol.json").read_text())
    rng=np.random.default_rng(cfg["seed"]); m=cfg["local_variants"]; n1=cfg["n1"]; n2=cfg["n2"]; reps=cfg["replicates_per_scenario"]
    native=[]
    for a in [.2,.6,.9,.97]:
        for h in [.0001,.001,.01]:
            s,other,err,_=native_fixture(a,h,.2)
            native.append({"ld_ar1":a,"local_h2":h,**{f"absolute_error_{k}":e for k,e in err.items()}})
            if err["p"]>1e-10 or err["var"]>1e-15 or err["rho"]>1e-12: raise ValueError("Native/reconstruction discordance")
    write_tsv("upstream_concordance.tsv",native)
    rows=[]
    for a in cfg["ar1_ld"]:
        r,d,v,l2=exact_ld(m,a); meanld=float(np.mean(l2))
        for h in cfg["local_h2_both_traits"]:
            av=d+n1*h*d*d/m; bv=d+n2*h*d*d/m
            for c in cfg["true_overlap_intercept"]:
                cross=c*d
                if np.any(av*bv-cross*cross<=0): raise ValueError("invalid generating covariance")
                analytic_var=float(np.sum(av*bv+cross*cross)/(n1*n2*meanld*meanld))
                vals={k:[] for k in ["rho","var","p","corr","oracle_moment_p","oracle_gls_p","oracle_gls_rho","selected_modes_for_variance"]}
                for start in range(0,reps,cfg["batch_size"]):
                    b=min(cfg["batch_size"],reps-start)
                    tx=np.sqrt(av)[None,None,:]*rng.standard_normal((b,cfg["context_blocks"],m))
                    ty=cross/av*tx+np.sqrt(bv-cross*cross/av)[None,None,:]*rng.standard_normal((b,cfg["context_blocks"],m))
                    hc1=(np.mean(tx*tx,axis=(1,2))-1)*m/(n1*meanld)
                    hc2=(np.mean(ty*ty,axis=(1,2))-1)*m/(n2*meanld)
                    est=reconstructed(tx[:,0,:],ty[:,0,:],d,n1,n2,hc1,hc2,meanld,c)
                    y=tx[:,0,:]*ty[:,0,:]-cross
                    vary=av*bv+cross*cross
                    den=np.sum(d**4/vary)
                    gls=m/math.sqrt(n1*n2)*np.sum(d*d*y/vary,axis=1)/den
                    glsvar=m*m/(n1*n2*den)
                    for k in vals:
                        value=(2*norm.sf(np.abs(est["rho"])/math.sqrt(analytic_var)) if k=="oracle_moment_p" else 2*norm.sf(np.abs(gls)/math.sqrt(glsvar)) if k=="oracle_gls_p" else gls if k=="oracle_gls_rho" else est[k])
                        vals[k].extend(value.tolist())
                out={k:np.asarray(x) for k,x in vals.items()}; valid=np.isfinite(out["p"])&(out["var"]>0)
                record={"scenario_id":f"ar{a}_h{h}_c{c}","classification":"TRUTH_KNOWN_SIMULATION","ld_ar1":a,"local_h2_both":h,"true_covariance":0,"true_intercept":c,"replicates":reps,"valid_estimates":int(valid.sum()),"invalid_estimates":int((~valid).sum()),"mean_estimated_covariance":float(np.mean(out["rho"])),"mean_covariance_mcse":float(np.std(out["rho"],ddof=1)/math.sqrt(reps)),"empirical_covariance_variance":float(np.var(out["rho"],ddof=1)),"mean_reported_variance":float(np.mean(out["var"][valid])),"exact_moment_variance":analytic_var,"empirical_over_reported_variance":float(np.var(out["rho"],ddof=1)/np.mean(out["var"][valid])),"exact_over_reported_variance":float(analytic_var/np.mean(out["var"][valid])),"numeric_correlations":int(np.isfinite(out["corr"]).sum()),"out_of_bounds_correlations":int((np.isfinite(out["corr"])&(np.abs(out["corr"])>1)).sum()),"median_selected_modes":float(np.median(out["selected_modes_for_variance"]))}
                for label in ["p","oracle_moment_p","oracle_gls_p"]:
                    for alpha in cfg["diagnostic_thresholds"]:
                        k=int((out[label]<=alpha).sum()); low,high=wilson(k,reps); prefix=f"{label}_alpha{alpha}"
                        record.update({prefix+"_rejects":k,prefix+"_rate":k/reps,prefix+"_wilson_lower":low,prefix+"_wilson_upper":high})
                rows.append(record); write_tsv("calibration_simulations.tsv",rows)
                print(record["scenario_id"],"variance_ratio",round(record["exact_over_reported_variance"],4),"nominal_rejections",record["p_alpha0.05_rejects"],flush=True)
    result={"classification":cfg["classification"],"scenarios_completed":len(rows),"total_replicates":len(rows)*reps,"native_function_concordance_fixtures":len(native),"upstream_source_sha256":SOURCE_SHA,"freeze_sha256":sha(HERE/"calibration_freeze.json"),"settings":cfg,"finished_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (HERE/"calibration_summary.json").write_text(json.dumps(result,indent=2)+"\n")

if __name__ == "__main__":
    import sys
    if sys.argv[1:]==["freeze"]: freeze()
    elif sys.argv[1:]==["run"]: run()
    else: raise SystemExit("Usage: calibration.py freeze|run")
