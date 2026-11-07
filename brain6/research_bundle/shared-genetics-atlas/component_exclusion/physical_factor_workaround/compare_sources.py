import json,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
reg=json.loads((p/'factor_feasibility_insomnia.json').read_text()); q=json.loads((p/'factor_feasibility_questionnaire.json').read_text())
assert reg['gate'].startswith('factor admissible') and q['gate'].startswith('factor admissible')
def rg_delete(a,b):
 base=p/(a+'_'+b+a+'.common.sumstats.gz_'+b+'.common.sumstats.gz')
 h1=np.loadtxt(str(base)+'.hsq1.delete');h2=np.loadtxt(str(base)+'.hsq2.delete');g=np.loadtxt(str(base)+'.gencov.delete')
 assert len(h1)==len(h2)==len(g)==200 and (h1>0).all() and (h2>0).all()
 return g/np.sqrt(h1*h2)
wg=rg_delete('walking','grip');wa=rg_delete('walking','activity');ga=rg_delete('grip','activity')
assert (wg>0).all() and (wa>0).all() and (ga>0).all()
lam=np.column_stack([np.sqrt(wg*wa/ga),np.sqrt(wg*ga/wa),np.sqrt(wa*ga/wg)])
block={}
for sleep,result in [('insomnia',reg),('questionnaire',q)]:
 rs=np.column_stack([rg_delete(t,sleep) for t in ('walking','grip','activity')]); weights=np.asarray(result['gls_weights']);block[sleep]=(rs/lam)@weights
delta=float(q['factor_insomnia_rg']-reg['factor_insomnia_rg']); d=block['questionnaire']-block['insomnia']; se=float(np.sqrt(199/200*np.sum((d-d.mean())**2)))
out={'questionnaire_minus_registry':delta,'paired_fixed_weight_block_se':se,'normal_p':math.erfc(abs(delta/se)/math.sqrt(2)),
     'unadjusted_95_ci':[delta-1.95996398454*se,delta+1.95996398454*se],
     'caveat':'Post-hoc comparison; weights fixed at full-data estimates. Does not include uncertainty from weight estimation or phenotype/cohort/source assumptions; not an independent-cohort comparison.'}
(p/'source_contrast.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
