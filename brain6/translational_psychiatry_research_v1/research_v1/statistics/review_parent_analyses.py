"""Separate post-result numerical review; never mutates parent frozen analysis files."""
import ast, hashlib, itertools, json, math, pathlib
from statistics import NormalDist
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from calibration import HERE, ROOT, sha, write_tsv

BASE=ROOT/'brain6/translational_psychiatry_research_v1'
DATA=BASE/'research_v1'

def exclusive_pair_logsum(l1,l2):
    # Sum all i!=j configurations in bounded blocks; no prefix/suffix subtraction.
    pieces=[]
    for start in range(0,len(l1),128):
        stop=min(len(l1),start+128)
        block=l1[start:stop,None]+l2[None,:]
        block[np.arange(stop-start),np.arange(start,stop)]=-np.inf
        pieces.append(logsumexp(block))
    return logsumexp(pieces)

def separate_posterior(l1,l2,p1,p2,p12):
    terms=np.array([0,math.log(p1)+logsumexp(l1),math.log(p2)+logsumexp(l2),math.log(p1)+math.log(p2)+exclusive_pair_logsum(l1,l2),math.log(p12)+logsumexp(l1+l2)])
    return np.exp(terms-logsumexp(terms))

def separate_logabf(beta,se,sd):
    variance=se*se; prior=sd*sd
    return .5*(np.log(variance)-np.log(variance+prior)+prior/(variance+prior)*beta*beta/variance)

def main():
    globalcfg=json.loads((DATA/'global_contrasts/protocol.json').read_text())
    globalinput=ROOT/'brain6/paper/final_package_v1/BRAIN6_FINAL_GLOBAL.tsv'
    assert sha(globalinput)==globalcfg['input_sha256']
    assert sha(BASE/'scripts/global_contrasts.py')==globalcfg['code_sha256']
    inherited=pd.read_csv(globalinput,sep='\t'); observed=pd.read_csv(DATA/'global_contrasts/contrasts.tsv',sep='\t')
    lookup={(r.sleep_trait,r.brain_disorder):r for r in inherited.itertuples()}; checks=[]
    critical=NormalDist().inv_cdf(1-.05/144); passing=0; eligible=0
    for r in observed.itertuples():
        x,y=lookup[r.sleep_trait,r.disorder1],lookup[r.sleep_trait,r.disorder2]
        delta=x.rg-y.rg; se=x.se+y.se
        marg=[]
        for v in [x,y]:
            marg.append(all(math.isfinite(getattr(v,k)) for k in ['rg','se','sleep_h2','sleep_h2_se','h2_disorder','h2_disorder_se']) and abs(v.rg)<=1 and min(v.se,v.sleep_h2,v.sleep_h2_se,v.h2_disorder,v.h2_disorder_se)>0 and v.sleep_h2/v.sleep_h2_se>=4 and v.h2_disorder/v.h2_disorder_se>=4)
        admitted=all(marg); eligible+=admitted
        # stdlib erfc independently computes the same normal tail; no parent imports.
        p=math.erfc(abs(delta)/se/math.sqrt(2)); rejected=admitted and p<=.05/72; passing+=rejected
        checks.append({'sleep_trait':r.sleep_trait,'disorder1':r.disorder1,'disorder2':r.disorder2,'reviewed_admitted':admitted,'reviewed_rejected':rejected,'delta_abs_error':abs(delta-r.difference),'SE_abs_error':abs(se-r.maximum_sampling_SE),'P_abs_error':abs(p-r.maximum_normal_model_P),'adjusted_P_abs_error':abs(min(1,72*p)-r.bonferroni72_P),'CI_lower_abs_error':abs(delta-critical*se-r.simultaneous_lower),'CI_upper_abs_error':abs(delta+critical*se-r.simultaneous_upper),'classification':'POST_RESULT_SEPARATE_IMPLEMENTATION_REVIEW'})
    expected={(t,a,b) for t in inherited.sleep_trait.unique() for a,b in itertools.combinations(globalcfg['disorders'],2)}
    assert len(checks)==72 and {(r['sleep_trait'],r['disorder1'],r['disorder2']) for r in checks}==expected
    assert eligible==72 and passing==16
    assert max(r[k] for r in checks for k in ['delta_abs_error','SE_abs_error','P_abs_error','adjusted_P_abs_error','CI_lower_abs_error','CI_upper_abs_error'])<1e-12
    write_tsv('parent_global_contrast_checks.tsv',checks)
    molcfg=json.loads((DATA/'molecular/sensitivity_protocol.json').read_text())
    amendment=json.loads((DATA/'molecular/serialization_amendment.json').read_text())
    assert sha(DATA/'molecular/frozen_analysis_code_v1.py.txt')==molcfg['code_sha256']==amendment['original_code_sha256']
    assert sha(BASE/'scripts/molecular_sensitivity.py')==amendment['corrected_code_sha256']
    assert sha(DATA/'molecular/coloc_5.2.3_source.R')==molcfg['native_reference_source_sha256']
    original=ast.parse((DATA/'molecular/frozen_analysis_code_v1.py.txt').read_text())
    updated=ast.parse((BASE/'scripts/molecular_sensitivity.py').read_text())
    unchanged_math=[]
    for name in ['estimate_sdy','logabf','posterior']:
        a=next(n for n in original.body if isinstance(n,ast.FunctionDef) and n.name==name)
        b=next(n for n in updated.body if isinstance(n,ast.FunctionDef) and n.name==name)
        assert ast.dump(a)==ast.dump(b); unchanged_math.append(name)
    mol=pd.read_csv(DATA/'molecular/sensitivity_results.tsv',sep='\t'); old=pd.read_csv(ROOT/'brain6/results/brain6_exploratory_functional_v1/coloc_results.tsv',sep='\t'); checks=[]; default=[]
    for p,digest in molcfg['files'].items():
        assert sha(ROOT/p)==digest
        d=pd.read_csv(ROOT/p,sep='\t')
        # Independent least-squares solver, rather than the parent's dot-products.
        predictor=(1/d.qtl_se.to_numpy()**2)[:,None]
        response=2*d.qtl_n.to_numpy()*d.qtl_maf.to_numpy()*(1-d.qtl_maf.to_numpy())
        sdy=float(math.sqrt(np.linalg.lstsq(predictor,response,rcond=None)[0][0]))
        l1=separate_logabf(d.gwas_beta.to_numpy(),d.gwas_se.to_numpy(),.2)
        for scale in molcfg['qtl_sdY_scale']:
            l2=separate_logabf(d.qtl_beta.to_numpy(),d.qtl_se.to_numpy(),.15*sdy*scale)
            # H3 pair-grid calculation cached across the five p12 values.
            distinct=exclusive_pair_logsum(l1,l2)
            base=np.array([0,math.log(molcfg['p1'])+logsumexp(l1),math.log(molcfg['p2'])+logsumexp(l2),math.log(molcfg['p1'])+math.log(molcfg['p2'])+distinct,logsumexp(l1+l2)])
            for p12 in molcfg['p12']:
                terms=base.copy();terms[4]+=math.log(p12); pp=np.exp(terms-logsumexp(terms))
                r=mol[(mol.input_path==p)&(mol.sdY_scale==scale)&(mol.p12==p12)].iloc[0]
                checks.append({'input_path':p,'sdY_scale':scale,'p12':p12,'reviewed_sdY':sdy,'sdY_abs_error':abs(sdy-r.estimated_sdY),'max_pp_abs_error':max(abs(pp[i]-r[f'pp_h{i}']) for i in range(5)),'reviewed_H4':pp[4],'H4_ge08':bool(pp[4]>=.8),'classification':'POST_RESULT_FULL_CONFIGURATION_SUM_REVIEW'})
                if scale==1 and p12==1e-5:
                    historical=old[old.input_path==p].iloc[0]
                    default.append({'input_path':p,'historical_default_max_pp_abs_error':max(abs(pp[i]-historical[f'pp_h{i}']) for i in range(5))})
    assert len(checks)==60 and not any(r['H4_ge08'] for r in checks)
    assert max(r['max_pp_abs_error'] for r in checks)<1e-10
    assert max(r['historical_default_max_pp_abs_error'] for r in default)<1e-10
    write_tsv('parent_molecular_configuration_checks.tsv',checks)
    write_tsv('parent_molecular_default_replay_checks.tsv',default)
    review=json.loads((DATA/'independent_statistics/review_settings.json').read_text())
    assert sha(BASE/'scripts/independent_moment_check.py')==review['code_sha256']
    existing=json.loads((DATA/'independent_statistics/independent_variance_results.json').read_text()); qc=[]
    m=review['m'];n1=review['n1'];n2=review['n2']
    for r in existing:
        ar,h,c=r['ld_ar1'],r['local_h2'],r['overlap_intercept']
        R=ar**np.abs(np.arange(m)[:,None]-np.arange(m)[None,:]);rr=R@R
        A=R+n1*h/m*rr;B=R+n2*h/m*rr;C=c*R
        Omega=np.block([[A,C],[C.T,B]])
        H=np.block([[np.zeros((m,m)),np.eye(m)/2],[np.eye(m)/2,np.zeros((m,m))]])
        denominator=n1*n2*(np.sum(R*R)/m)**2
        dense_quad=2*np.trace((H@Omega)@(H@Omega))/denominator
        relative=abs(dense_quad-r['trace_based_variance'])/dense_quad
        qc.append({'ld_ar1':ar,'local_h2':h,'intercept':c,'separate_dense_quadratic_form_variance':dense_quad,'parent_trace_variance':r['trace_based_variance'],'relative_error':relative,'parent_MC_variance':r['empirical_variance'],'parent_reported_MCSE':r['variance_mcse'],'analytic_vs_parent_MC_in_MCSE':abs(dense_quad-r['empirical_variance'])/r['variance_mcse'],'classification':'POST_RESULT_DENSE_QUADRATIC_FORM_REVIEW'})
        assert relative<1e-12
    write_tsv('parent_independent_variance_checks.tsv',qc)
    global_review=pd.read_csv(HERE/'parent_global_contrast_checks.tsv',sep='\t')
    global_error_columns=[c for c in global_review if c.endswith('abs_error')]
    summary={'classification':'SEPARATE_POST_RESULT_AUTOMATED_ADVERSARIAL_REVIEW','global_contrasts_reviewed':72,'global_unique_marginal_rows':48,'global_admitted':eligible,'global_conservative_differences':passing,'max_global_numeric_error':float(global_review[global_error_columns].to_numpy().max()),'molecular_conditions_reviewed':60,'max_molecular_pp_error':float(max(r['max_pp_abs_error'] for r in checks)),'max_molecular_H4':float(max(r['reviewed_H4'] for r in checks)),'molecular_H4_ge08':0,'unchanged_molecular_math_AST_functions':unchanged_math,'independent_variance_cases':2,'max_dense_quadratic_form_relative_error':max(r['relative_error'] for r in qc),'max_analytic_vs_parent_MC_in_MCSE':max(r['analytic_vs_parent_MC_in_MCSE'] for r in qc),'not_independent_human_validation':True,'not_independent_cohort_replication':True}
    (HERE/'parent_analysis_review_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
