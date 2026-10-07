"""Post-result deterministic Gaussian-product tail checks; does not modify frozen experiments."""
import json
import numpy as np
from calibration import HERE, exact_ld, write_tsv

def product_cumulants(a,b,c):
    # log M_XY(t)=-1/2 log[1-2ct-(ab-c^2)t^2]; subtracting mean leaves k2/k4 unchanged.
    return a*b+c*c, 6*(a*a*b*b+6*a*b*c*c+c**4)

def main():
    cfg=json.loads((HERE/'calibration_protocol.json').read_text()); m=cfg['local_variants']; rows=[]
    for ar in cfg['ar1_ld']:
        _,d,_,l2=exact_ld(m,ar)
        for h in cfg['local_h2_both_traits']:
            a=d+cfg['n1']*h*d*d/m; b=d+cfg['n2']*h*d*d/m
            for overlap in cfg['true_overlap_intercept']:
                v,k4=product_cumulants(a,b,overlap*d)
                w=d*d/v; glsv=w*w*v
                h1se=float(np.sqrt(2*np.sum(a*a))/(cfg['n1']*np.mean(l2)))
                h2se=float(np.sqrt(2*np.sum(b*b))/(cfg['n2']*np.mean(l2)))
                rows.append({'scenario_id':f'ar{ar}_h{h}_c{overlap}','classification':'POST_RESULT_DETERMINISTIC_MODEL_TAIL_DIAGNOSTIC','moment_variance_effective_mode_count':float(v.sum()**2/(v*v).sum()),'moment_largest_mode_variance_share':float(v.max()/v.sum()),'moment_exact_excess_kurtosis':float(k4.sum()/v.sum()**2),'oracle_gls_variance_effective_mode_count':float(glsv.sum()**2/(glsv*glsv).sum()),'oracle_gls_largest_mode_variance_share':float(glsv.max()/glsv.sum()),'oracle_gls_exact_excess_kurtosis':float((w**4*k4).sum()/glsv.sum()**2),'true_local_h2_both_traits':h,'moment_local_h2_1_exact_se':h1se,'moment_local_h2_2_exact_se':h2se,'local_h2_1_relative_se':h1se/h,'local_h2_2_relative_se':h2se/h,'caveat':'Oracle normal Wald P is a benchmark, not exact finite-sample calibration; positive kurtosis/heavy tails remain. Synthetic h2 uncertainty is not an empirical Brain6 confidence interval.'})
    write_tsv('gaussian_product_tail_diagnostics.tsv',rows)

if __name__=='__main__': main()
