"""Independent single-signal ABF diagnostic on four inherited molecular slices."""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'brain6/results/brain6_exploratory_functional_v1'
OUT=ROOT/'brain6/translational_psychiatry_research_v1/research_v1/molecular'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def estimate_sdy(v,maf,n):
    x=1/v; y=2*n*maf*(1-maf)
    return np.sqrt(np.dot(x,y)/np.dot(x,x))
def logabf(beta,se,prior_sd):
    v=se**2; r=prior_sd**2/(v+prior_sd**2)
    return .5*(np.log1p(-r)+r*(beta/se)**2)
def posterior(l1,l2,p1,p2,p12):
    # Independently sum all distinct-variant products without subtracting near-equal exponentials.
    a=np.exp(l1-np.max(l1)); b=np.exp(l2-np.max(l2))
    prefix=np.cumsum(b)-b; suffix=np.cumsum(b[::-1])[::-1]-b
    distinct=np.dot(a,prefix+suffix)
    logs=np.array([0,np.log(p1)+logsumexp(l1),np.log(p2)+logsumexp(l2),np.log(p1*p2)+np.log(distinct)+np.max(l1)+np.max(l2),np.log(p12)+logsumexp(l1+l2)])
    return np.exp(logs-logsumexp(logs))
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=sorted((OLD/'coloc_inputs').glob('*.gz')); assert len(paths)==4
    cfgpath=OUT/'sensitivity_protocol.json'
    if sys.argv[1:]==['freeze']:
        if cfgpath.exists():raise ValueError('Do not overwrite freeze')
        cfg={'id':'MOLECULAR_ABF_DIAGNOSTIC_V1','classification':'HISTORICAL_INPUT_MODEL_SENSITIVITY_NOT_NEW_MECHANISM','files':{str(p.relative_to(ROOT)):sha(p) for p in paths},'code_sha256':sha(Path(__file__)),'native_reference_source_sha256':sha(OUT/'coloc_5.2.3_source.R'),'p1':1e-4,'p2':1e-4,'p12':[1e-6,5e-6,1e-5,5e-5,1e-4],'qtl_sdY_scale':[.5,1.,2.],'family':'All four originally eligible ADHD-component molecular contexts;60posterior conditions','default_replay_tolerance':1e-10,'decision':'No causal/gene/pair promotion; all prior conditions retained','limitations':['Single causal variant approximation','Same discovery ADHD GWAS','No insomnia component molecular test','sdY estimated under homoscedastic regression assumption','Source molecular phenotype scale not independently validated','Alleles already lifted/matched upstream; cannot independently recover omitted strands from these slices']}
        cfgpath.write_text(json.dumps(cfg,indent=2)+'\n');return
    if sys.argv[1:]!=['run']:raise SystemExit('freeze|run')
    cfg=json.loads(cfgpath.read_text())
    if cfg['code_sha256']!=sha(Path(__file__)):
        amendment=json.loads((OUT/'serialization_amendment.json').read_text())
        assert amendment['original_code_sha256']==cfg['code_sha256']
        assert amendment['corrected_code_sha256']==sha(Path(__file__))
    expected=pd.read_csv(OLD/'coloc_results.tsv',sep='\t'); results=[];checks=[]
    for path in paths:
        assert cfg['files'][str(path.relative_to(ROOT))]==sha(path)
        d=pd.read_csv(path,sep='\t'); assert not d.variant.duplicated().any()
        assert np.all(np.isfinite(d.select_dtypes('number'))) and np.all(d.gwas_se>0) and np.all(d.qtl_se>0)
        assert np.all((d.qtl_maf>0)&(d.qtl_maf<=.5)) and np.all(d.qtl_n>0)
        old=expected[expected.input_path==str(path.relative_to(ROOT))].iloc[0]
        sd=estimate_sdy(d.qtl_se.to_numpy()**2,d.qtl_maf.to_numpy(),d.qtl_n.to_numpy())
        l1=logabf(d.gwas_beta.to_numpy(),d.gwas_se.to_numpy(),.2)
        for scale in cfg['qtl_sdY_scale']:
            l2=logabf(d.qtl_beta.to_numpy(),d.qtl_se.to_numpy(),.15*sd*scale)
            for p12 in cfg['p12']:
                pp=posterior(l1,l2,cfg['p1'],cfg['p2'],p12)
                row={'input_path':str(path.relative_to(ROOT)),'gene_id':old.gene_id,'molecular_trait_id':old.molecular_trait_id,'context':old.qtl_context,'n_shared':len(d),'estimated_sdY':sd,'sdY_scale':scale,'p12':p12,**{f'pp_h{i}':p for i,p in enumerate(pp)},'status':cfg['classification'],'input_sha256':sha(path)}
                results.append(row)
                if scale==1 and p12==1e-5:
                    err=max(abs(pp[i]-float(old[f'pp_h{i}'])) for i in range(5))
                    checks.append({'input_path':str(path.relative_to(ROOT)),'maximum_posterior_absolute_error':err,'status':'PASS' if err<=cfg['default_replay_tolerance'] else 'FAIL'})
                    if err>cfg['default_replay_tolerance']:raise ValueError('Independent posterior mismatch')
    for name,rows in [('sensitivity_results.tsv',results),('independent_posterior_checks.tsv',checks)]:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    summary={'actual_molecular_inputs':4,'posterior_conditions':len(results),'default_replays':len(checks),'max_H4_across_all_conditions':float(max(r['pp_h4'] for r in results)),'conditions_H4_ge08':int(sum(r['pp_h4']>=.8 for r in results)),'classification':cfg['classification']}
    (OUT/'sensitivity_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
if __name__=='__main__':main()
