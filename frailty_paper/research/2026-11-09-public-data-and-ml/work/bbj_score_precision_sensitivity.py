import runpy,json,math,csv
import numpy as np
m=runpy.run_path('outputs/stage6_bbj_insomnia_matched14_z.py')
selected=m['selected'];bbj=m['bbj'];chrom=m['chrom'];x=m['x'];mask=m['mask']; scale=m['scale'];z=np.array([c['beta_over_se'] for c in m['components']]);se=np.array([float(bbj[r['variant']]['se']) for r in selected]);scaled_se=se*scale
R=np.corrcoef(x[mask],rowvar=False);R=np.where(np.equal.outer(chrom,chrom),R,0)
def calc(a):return float(a@z/math.sqrt(a@R@a))
frequency=calc(scale);precision=calc(1/se)
assert abs(frequency-m['result']['JPT_n104_beta_over_se']['z'])<1e-10
rows=[{'variant':r['variant'],'allele_frequency':float(c['freq']),'standard_error':float(e),'se_times_genotype_sd':float(e*s)} for r,c,e,s in zip(selected,m['components'],se,scale)]
with open('outputs/stage6_bbj_14_precision_components.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
r={'n_markers':len(selected),'se_times_genotype_sd_min':float(scaled_se.min()),'se_times_genotype_sd_median':float(np.median(scaled_se)),'se_times_genotype_sd_max':float(scaled_se.max()),'coefficient_of_variation':float(np.std(scaled_se)/np.mean(scaled_se)),'frequency_scaled_z':frequency,'inverse_se_scaled_z':precision,'method':'Score weights assumed equal; compare Z weighting by sqrt(2p(1-p)) versus 1/SE, with same target marginal beta/SE and JPT reference LD; cross-chromosome correlations zero.','caveat':'Postselected exploratory analysis. Inverse-SE assumes these marginal standard errors and external LD approximate the target score-statistic covariance; source imbalance, imputation, sample overlap, low case counts and target LD can invalidate either approximation. Not a patient-level score test, association confidence interval, or clinical validation.'}
open('outputs/stage6_bbj_14_precision_sensitivity.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
