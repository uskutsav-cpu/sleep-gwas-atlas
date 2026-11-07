import csv,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent;atlas=p.parents[1];raw=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run/predictions.tsv';one=atlas/'global_audit/mvp_adhd_frailty_partial/results.tsv';two=atlas/'global_audit/mvp_adhd_mdd_frailty_partial/results.tsv'
def read(f,key):
 with f.open(newline='') as h:return {r[key]:r for r in csv.DictReader(h,delimiter='\t')}
a=read(raw,'target');b=read(one,'frailty');c=read(two,'frailty');order=['FI','FFS','HFRS'];assert set(a)==set(b)==set(c)==set(order)
fig,axes=plt.subplots(1,3,figsize=(11.4,4.2),sharey=True,layout='constrained');panels=[('Genome-wide r_g',a,'observed_mvp_rg','observed_mvp_se','#325c81'),('Partial r_g given ADHD',b,'partial_genetic_rg','paired_jackknife_se','#784b83'),('Partial r_g given ADHD + depression',c,'partial_rg_given_adhd_mdd','paired_jackknife_se','#9c5b32')]
for ax,(title,data,point,err,color) in zip(axes,panels):
 for i,name in enumerate(order):
  x=float(data[name][point]);se=float(data[name][err]);ax.errorbar(x,2-i,xerr=1.96*se,fmt='o',color=color,ecolor=color,elinewidth=1.8,capsize=3,markersize=6,zorder=3)
 ax.axvline(0,color='#777777',lw=.8,ls='--');ax.set_xlim(-.15,.75);ax.set_ylim(-.65,2.65);ax.set_yticks([2,1,0],['Frailty index (FI)','Fried score (FFS)','Hospital score (HFRS)']);ax.set_xlabel('Genetic correlation');ax.set_title(title,fontsize=11,fontweight='bold');ax.grid(axis='x',color='#e8e8e8',lw=.6);ax.spines[['top','right']].set_visible(False)
fig.suptitle('Clinical insomnia and frailty: estimates on shared markers',fontsize=14,fontweight='bold');fig.text(.5,-.055,'Points and approximate 95% intervals; partial correlations are genetic projections, not causal or diagnosis-adjusted effects.\nThe scores contain sleep, fatigue, exhaustion or somnolence items; source overlap and model calibration remain uncertain.',ha='center',fontsize=8.3,color='#444444')
for fmt in ('png','svg'):fig.savefig(p/('figure.'+fmt),dpi=180,bbox_inches='tight',facecolor='white')
(p/'receipt.json').write_text(json.dumps({'source_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (raw,one,two)},'snps':991707,'interval':'point ± 1.96 reported LDSC or paired-jackknife SE; approximate and unadjusted','note':'post-hoc estimates, differing estimands; not individual adjustment or causal effects'},indent=2)+'\n');print('saved figure')
