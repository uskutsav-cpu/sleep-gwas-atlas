"""Post-result independent analytic/Monte Carlo falsification of a variance mismatch."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'brain6/translational_psychiatry_research_v1/research_v1/independent_statistics'
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # These fixed counterexamples were chosen after agent diagnostic results. No discovery claim.
    settings={'classification':'POST_RESULT_AUTOMATED_INDEPENDENT_IMPLEMENTATION_REVIEW','scenarios':[[.2,.0001,0],[.97,.01,.2]],'m':160,'n1':386533,'n2':225534,'replicates':100000,'seed':2026100751,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'formula':'trace(Sigma1 Sigma2)+trace(C C) for Gaussian z1 dot z2; compare source-specific analytic eigenvariance','not_empirical_GWAS_calibration':True}
    (OUT/'review_settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    rng=np.random.default_rng(settings['seed']);rows=[];m=settings['m'];n1=settings['n1'];n2=settings['n2']
    for ar,h,c in settings['scenarios']:
        R=ar**np.abs(np.arange(m)[:,None]-np.arange(m)[None,:]);D,U=np.linalg.eigh(R)
        Sigma1=R+n1*h/m*(R@R);Sigma2=R+n2*h/m*(R@R);C=c*R
        denominator=n1*n2*(np.trace(R@R)/m)**2
        var_trace=(np.trace(Sigma1@Sigma2)+np.trace(C@C))/denominator
        a=D+n1*h/m*D**2;b=D+n2*h/m*D**2;cc=c*D
        var_modes=np.sum(a*b+cc**2)/denominator
        moments=[]
        for start in range(0,settings['replicates'],500):
            x=np.sqrt(a)*rng.normal(size=(500,m));y=cc/a*x+np.sqrt(b-cc**2/a)*rng.normal(size=(500,m))
            moments.extend(((x*y).sum(axis=1)-c*m)/(np.sqrt(n1*n2)*np.trace(R@R)/m))
        values=np.array(moments); empirical=np.var(values,ddof=1);mcse=np.std(values**2,ddof=1)/np.sqrt(len(values))
        rows.append({'ld_ar1':ar,'local_h2':h,'overlap_intercept':c,'trace_based_variance':var_trace,'eigen_based_variance':var_modes,'relative_analytic_error':abs(var_trace-var_modes)/var_trace,'empirical_variance':empirical,'variance_mcse':mcse,'analytic_vs_MC_deviations':abs(empirical-var_trace)/mcse,'replicates':len(values)})
        assert abs(var_trace-var_modes)/var_trace<1e-12
    (OUT/'independent_variance_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows)
if __name__=='__main__':main()
