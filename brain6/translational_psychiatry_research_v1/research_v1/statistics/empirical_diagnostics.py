"""Retrospective source-specific diagnostics; never recalculate empirical genetic results from synthetic LD."""
import datetime, hashlib, json, math, pathlib, shutil, subprocess, tempfile
import numpy as np
import pandas as pd
from scipy.stats import norm
from calibration import HERE, ROOT, native_fixture, reconstructed, write_tsv, sha

def main():
    raw=ROOT/'brain6/results/brain6_alternative_local_validation_v1/genome_wide_blocks.tsv'
    source=ROOT/'brain6/results/brain6_alternative_local_validation_v1/protocol_freeze.json'
    cfg=json.loads(source.read_text()); x=pd.read_csv(raw,sep='\t'); x=x[x.analysis_status=='ESTIMATED'].copy()
    critical=float(norm.isf(cfg['alpha_block']/2))
    rows=[]
    for r in x.to_dict('records'):
        positive=r['h2_1']>0 and r['h2_2']>0
        bound=math.sqrt(r['h2_1']*r['h2_2']) if positive else float('nan')
        se=math.sqrt(r['var']); z=abs(r['rho'])/se
        rows.append({'pair_id':r['pair_id'],'chr':int(r['chr']),'start':int(r['start']),'end':int(r['end']),'m':int(r['m']),'rho':r['rho'],'h2_1':r['h2_1'],'h2_2':r['h2_2'],'corr':r['corr'],'covariance_se':se,'both_h2_positive':positive,'derived_correlation_out_of_bounds':bool(positive and abs(r['corr'])>1),'moment_parameter_matrix_determinant':r['h2_1']*r['h2_2']-r['rho']*r['rho'],'physical_bound_using_point_h2':bound,'physical_excess_over_covariance_se_DIAGNOSTIC_NOT_TEST':(abs(r['rho'])-bound)/se if positive else float('nan'),'source_wald_abs_z':z,'source_family_significant':bool(r['p']<=cfg['alpha_block']),'minimum_se_multiplier_to_lose_significance':z/critical if r['p']<=cfg['alpha_block'] else float('nan'),'fixed_source_wald_family_ci_lower':r['rho']-critical*se,'fixed_source_wald_family_ci_upper':r['rho']+critical*se,'classification':'RETROSPECTIVE_SOURCE_ARITHMETIC_UNCERTIFIED_MODEL','source_sha256':sha(raw)})
    write_tsv('empirical_block_diagnostics.tsv',rows)
    s=[r for r in rows if r['source_family_significant']]; write_tsv('significant_block_fragility.tsv',s)
    pair=[]
    for p,g in pd.DataFrame(rows).groupby('pair_id'):
        pair.append({'pair_id':p,'estimated':len(g),'both_h2_positive':int(g.both_h2_positive.sum()),'negative_h2_1':int((g.h2_1<0).sum()),'negative_h2_2':int((g.h2_2<0).sum()),'corr_missing':int(g['corr'].isna().sum()),'corr_out_of_bounds':int(g.derived_correlation_out_of_bounds.sum()),'source_family_significant':int(g.source_family_significant.sum()),'significant_corr_out_of_bounds':int((g.source_family_significant&g.derived_correlation_out_of_bounds).sum())})
    write_tsv('empirical_pair_diagnostics.tsv',pair)
    # Identity checks use the exact archived numerical routine on synthetic fixtures.
    scale=[]
    baseline=native_fixture()[0].iloc[0]
    for s1,s2 in [(.5,.5),(.5,2.),(2.,.5),(2.,2.)]:
        n1=int(386533*s1); n2=int(225534*s2)
        # Hold synthetic signed Z and LD fixed, rather than regenerate changed GWAS.
        _,_,_,(tx,ty,d,meanld,hc1,hc2)=native_fixture()
        est=reconstructed(tx[:1],ty[:1],d,n1,n2,np.array([hc1*386533/n1]),np.array([hc2*225534/n2]),meanld,.2)
        rho_factor=math.sqrt(386533*225534/(n1*n2)); var_factor=386533*225534/(n1*n2)
        scale.append({'n1_multiplier':n1/386533,'n2_multiplier':n2/225534,'rho_observed_ratio':est['rho'][0]/baseline.rho,'rho_expected_ratio':rho_factor,'variance_observed_ratio':est['var'][0]/baseline['var'],'variance_expected_ratio':var_factor,'p_absolute_difference':abs(est['p'][0]-baseline.p),'corr_absolute_difference':abs(est['corr'][0]-baseline['corr']),'classification':'SYNTHETIC_FIXED_Z_N_RESCALING_IDENTITY_NOT_SOURCE_N_VALIDATION'})
    write_tsv('scalar_n_invariance.tsv',scale)
    patch=ROOT/'brain6/results/brain6_alternative_local_validation_v1/SUPERGNOVA_strict_runtime.patch'
    expected={pathlib.Path(r['path']).name:r['sha256'] for r in cfg['frozen_files'] if '/SUPERGNOVA_source/' in r['path'] and '/._' not in r['path']}
    receipts=[]
    with tempfile.TemporaryDirectory(prefix='patched-method-',dir=HERE) as tmp:
        for name in ['calculate.py','prep.py','supergnova.py']: shutil.copy2(HERE/'source'/name,pathlib.Path(tmp)/name)
        p=subprocess.run(['patch','--batch','-p1','-i',str(patch)],cwd=tmp,capture_output=True,text=True)
        for path in sorted((HERE/'source').glob('*.py')):
            recreated=pathlib.Path(tmp)/path.name
            digest=sha(recreated if recreated.exists() else path)
            receipts.append({'file':path.name,'public_pinned_sha256':sha(path),'frozen_historical_sha256':expected.get(path.name,''),'reconstructed_historical_sha256':digest,'status':'EXACT_MATCH' if expected.get(path.name)==digest else 'MISMATCH_OR_UNLISTED','patch_exit_code':p.returncode,'patch_sha256':sha(patch)})
        patch_output=p.stdout+'\n'+p.stderr
    write_tsv('historical_source_identity.tsv',receipts)
    (HERE/'historical_patch_reconstruction.log').write_text(patch_output)
    summary={'classification':'RETROSPECTIVE_SOURCE_DIAGNOSTICS','source_sha256':sha(raw),'frozen_protocol_sha256':sha(source),'family_denominator':8465,'alpha':cfg['alpha_block'],'critical_abs_z':critical,'estimated_blocks':len(x),'summary_by_pair':pair,'historical_source_reconstruction':receipts,'intercept_and_local_h2_uncertainty':'Original meanLD, LD eigenvalues, genome-wide intercept/variance and joint local-h2 uncertainty absent. No empirical recalibration or determinant significance test performed.','N_identity_interpretation':'For fixed Z and reference and scalar N changes, N*h products are invariant; covariance scales inverse sqrt(N1*N2), variance inverse N1*N2, while Wald P and derived correlation are invariant. This does not validate correct N semantics or variant-varying sample contributions.','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (HERE/'empirical_diagnostics_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__': main()
