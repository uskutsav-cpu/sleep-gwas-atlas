#!/usr/bin/env python3
"""Deterministic vector figures, only from archived audited numerical tables."""
import os,json,csv,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/brain6-continuation-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
from audit_evidence import ROOT,AREA
M=AREA/'MANUSCRIPT';F=M/'figures';F.mkdir(exist_ok=True)
P=ROOT/'brain6/paper/final_package_v1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'brain6-20261007'})
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight',metadata={'Creator':'Brain6 deterministic table audit','CreationDate':None,'ModDate':None});fig.savefig(F/(name+'.svg'),bbox_inches='tight',metadata={'Date':None});plt.close(fig)
g=pd.read_csv(P/'BRAIN6_FINAL_GLOBAL.tsv',sep='\t');g=g.rename(columns={'brain_disorder':'disorder','significance_under_original_396_family':'significant_fdr'});sleep=list(dict.fromkeys(g.sleep_trait));dis=['adhd','mdd','scz','bipolar','parkinson','alzheimer'];print('disorders',g.disorder.unique());actual=set(g.disorder)
dis=[x for x in dis if x in actual]+sorted(actual-set(dis));pivot=g.pivot(index='sleep_trait',columns='disorder',values='rg').reindex(index=sleep,columns=dis)
fig,ax=plt.subplots(figsize=(7.2,6.3));v=max(.5,np.nanmax(np.abs(pivot.values)));im=ax.imshow(pivot.values,cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto');ax.set_xticks(range(len(dis)),['ADHD' if x=='adhd' else 'MDD' if x=='mdd' else 'SCZ' if x=='scz' else 'Bipolar' if x=='bipolar' else 'PD' if x=='parkinson' else 'AD' for x in dis]);ax.set_yticks(range(len(sleep)),[{'sleepdur':'Sleep duration','shortsleep':'Short sleep','longsleep':'Long sleep','sleepiness':'Daytime sleepiness','accel_sleep_duration':'Accelerometer sleep duration'}.get(x,x.replace('_',' ').capitalize()) for x in sleep]);
for i,st in enumerate(sleep):
 for j,di in enumerate(dis):
  row=g[(g.sleep_trait==st)&(g.disorder==di)].iloc[0]
  sig=str(row['significant_fdr']).lower() in ['true','1'];ax.text(j,i,f"{row.rg:.2f}"+('*' if sig else ''),ha='center',va='center',color='white' if abs(row.rg)>.65*v else 'black',fontsize=8)
ax.set_title('Global map: inherited estimates and original atlas FDR',pad=12);fig.colorbar(im,ax=ax,label='Recorded genetic correlation',shrink=.75);fig.text(.1,.01,'* Original 396-test FDR significance; selected display is not a new family.',fontsize=8);fig.tight_layout(rect=(0,.025,1,1));save(fig,'figure1_global_map');g.to_csv(F/'figure1_source.tsv',sep='\t',index=False)
b=json.loads((AREA/'qc/baseline_verification.json').read_text());counts=b['lava_not_run_by_trait'];order=sorted(counts,key=counts.get,reverse=True);fig,ax=plt.subplots(figsize=(7.2,4.0));ax.barh([{'longsleep':'Long sleep','insomnia':'Insomnia','parkinson':'PD','mdd':'MDD','adhd':'ADHD','bipolar':'Bipolar','scz':'SCZ'}[x] for x in order],[counts[x] for x in order],color='#366a8c');ax.invert_yaxis();ax.axvline(873,color='#b04444',ls='--',label='Whole-family ceiling: 873');
for i,x in enumerate(order):ax.text(counts[x]+12,i,str(counts[x]),va='center')
ax.set_xlim(0,1500);ax.set_xlabel('NOT_RUN cells (unestimated, not tested nulls)');ax.set_title('Canonical LAVA: 3,720 NOT_RUN / 17,465 planned');ax.legend(loc='lower right',frameon=False);fig.tight_layout();save(fig,'figure2_lava_qc');pd.DataFrame({'trait':order,'not_run':[counts[x] for x in order]}).to_csv(F/'figure2_source.tsv',sep='\t',index=False)
z=pd.read_csv(AREA/'qc/significant_secondary_blocks.tsv',sep='\t');fig,ax=plt.subplots(figsize=(7.2,3.8));labels=[r.pair_id.replace('insomnia__','Insomnia–').upper()+f'\nchr{r.chr}:{r.start/1e6:.2f}–{r.end/1e6:.2f} Mb' for _,r in z.iterrows()];x=z.rho.values*1e4;err=1.96*np.sqrt(z['var'].values)*1e4;ax.errorbar(x,np.arange(len(z)),xerr=err,fmt='o',color='#366a8c',capsize=4);ax.set_yticks(range(len(z)),labels);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('Observed-scale local covariance × 10⁴ (95% Wald interval)');ax.set_title('Three secondary family-significant blocks');fig.text(.03,.025,'Same discovery GWAS; corrected across all 8,465 slots.\nChr11 ADHD derived correlation = 1.1037 (outside bounds).',fontsize=8);fig.tight_layout(rect=(0,.12,1,1));save(fig,'figure3_secondary_covariance');z.to_csv(F/'figure3_source.tsv',sep='\t',index=False)
c=pd.read_csv(AREA/'qc/coloc_prior_replay.tsv',sep='\t');print('coloc cols',list(c));
# Use audited table field names, not posterior reconstruction inside plotting.
defcol=[k for k in c if 'default' in k.lower() and ('h4' in k.lower() or 'pp' in k.lower())][0];lowcol='low_prior_H4_replayed'
labels=[relabel.split('_chr')[-1].replace('_',':',1).replace('_','–') for relabel in c.candidate_locus_id];fig,ax=plt.subplots(figsize=(7.2,4.4));yy=np.arange(len(c));ax.hlines(yy,c[lowcol],c[defcol],color='#999999');ax.scatter(c[defcol],yy,label='Default p12=10⁻⁵',color='#366a8c');ax.scatter(c[lowcol],yy,label='Lower p12=10⁻⁶',color='#b04444',marker='s');ax.set_yticks(yy,['chr'+s for s in labels]);ax.invert_yaxis();ax.axvline(.8,color='gray',ls='--',lw=1);ax.set_xlim(0,1.05);ax.set_xlabel('PP.H4 under single-causal-variant model');ax.set_title('Exploratory prior sensitivity: all six estimated intervals');ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.19),ncol=2);fig.tight_layout(rect=(0,.07,1,1));save(fig,'figure4_coloc_priors');c.to_csv(F/'figure4_source.tsv',sep='\t',index=False)
manifest=[]
for p in sorted(F.glob('*')):manifest.append({'file':str(p.relative_to(AREA)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(AREA/'qc/figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
