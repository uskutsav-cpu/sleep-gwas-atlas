#!/usr/bin/env python3
"""Exact stock display audit; formats a few hundred numbers, runs no estimator."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SPEC=importlib.util.spec_from_file_location('independent_subset_util',Path(__file__).with_name('independent_core_subset_checker_v4_3.py'))
U=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(U)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    plan=U.source_json(P/'manifests/ssd_native_execution_plan_v4_3.json')
    native=Path(plan['ssd_support_package'])/'native/core_reproduction_v1'
    capture=U.source_json(native/'core_rg_insomnia.full_precision.json')
    tokens=U.correlation_tokens(native/'core_rg_insomnia.log')
    fmt_inputs=[e['p'] for e in capture['estimates']]
    fields={'rg':('rg_ratio',None),'se':('rg_se',None),'z':('z',None),'p':('p',None),'h2_obs':('tot','hsq2'),'h2_obs_se':('tot_se','hsq2'),'h2_int':('intercept','hsq2'),'h2_int_se':('intercept_se','hsq2'),'gcov_int':('intercept','gencov'),'gcov_int_se':('intercept_se','gencov')}
    h2caps=[]
    for ident in ['insomnia','ms','melanoma','ldl','sleepdur','t2d','parkinson']:
        c=U.source_json(native/('core_h2_'+ident+'.full_precision.json'));e=c['estimates'][0]
        factor=U.scale_factor(c['arguments']['samp_prev'],c['arguments']['pop_prev'])
        q=[e['tot']*factor,e['tot_se']*factor,e['intercept'],e['intercept_se']]
        fmt_inputs+=q;h2caps.append((ident,c,q,U.h2_tokens(native/('core_h2_'+ident+'.log'))))
    # The formatter implements stock matrix scalar display with the pinned NumPy,
    # but does not import LDSC, Pandas/SciPy, source GWAS or regression functions.
    code='import json,sys,numpy as np;np.set_printoptions(linewidth=1000,precision=4);x=json.load(sys.stdin);print(json.dumps(dict(numpy=np.__version__,values=[str(np.matrix(v)).replace("[","").replace("]","").strip() for v in x])))'
    python=next(k for k in plan['dependencies_sha256'] if k.endswith('/.ldsc-env/bin/python'))
    assert sha(Path(python))==plan['dependencies_sha256'][python]
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','TMPDIR':str(Path(plan['ssd_support_package']).parents[1]/'tmp')}
    result=subprocess.run([python,'-B','-c',code],input=json.dumps(fmt_inputs),capture_output=True,text=True,check=True,env=env)
    formatted=json.loads(result.stdout);assert formatted['numpy']=='1.21.5'
    pairs=[]
    for i,e in enumerate(capture['estimates']):
        key=(U.trait(e['p1']),U.trait(e['p2']));n=tokens[key];checks={}
        for display,(field,component) in fields.items():
            v=e[component][field] if component else e[field]
            checks[display]=format(v,'.4f')==n[display]
        checks['scalar_numpy_P']=formatted['values'][i]==n['scalar_P']
        pairs.append({'pair_id':'__'.join(key),'exact_display_checks':checks,'all_pass':all(checks.values()),'scalar_P_native':n['scalar_P'],'scalar_P_independently_formatted':formatted['values'][i]})
    h2=[]
    for i,(ident,c,q,n) in enumerate(h2caps):
        outputs=formatted['values'][33+4*i:33+4*(i+1)]
        checks={key:actual==n[key] for key,actual in zip(['h2','se','intercept','intercept_se'],outputs)}
        h2.append({'trait':ident,'exact_display_checks':checks,'all_pass':all(checks.values()),'independently_formatted':outputs})
    assert all(sha(Path(s))==m['sha256'] for s,m in U.SEEN.items())
    receipt={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Already completed insomnia33 and selected7h2 jobs; formatting-only Python subprocess, no estimator','candidate_collator_imported_or_executed':False,'exact_summary_4f_checks':330,'exact_scalar_numpy_P_checks':33,'exact_standalone_h2_scalar_checks':28,'all_checks_pass':all(r['all_pass'] for r in pairs+h2),'numpy_formatter':formatted['numpy'],'numpy_formatter_code':code,'formatter_python_sha256':sha(Path(python)),'formatter_elapsed_resource_scope':'61single scalar matrix displays; no array/data fitting','pairs':pairs,'h2':h2,'consumed_small_input_hashes':U.SEEN,'checker_sha256':sha(Path(__file__)),'subset_util_sha256':sha(Path(__file__).with_name('independent_core_subset_checker_v4_3.py')),'cross_estimator_block_alignment_certified':False}
    target=P/'reviews/independent_native_display_receipt_v4.json';assert not target.exists();target.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['pairs','h2','consumed_small_input_hashes','numpy_formatter_code']},indent=2))

if __name__=='__main__':main()
