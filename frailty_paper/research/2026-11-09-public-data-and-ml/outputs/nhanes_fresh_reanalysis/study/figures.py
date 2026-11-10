#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;res=pd.read_csv(P/'supplement_survey_models.csv');des=pd.read_csv(P/'supplement_weighted_descriptives.csv');d=pd.read_csv(P/'finalized_participant_data.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.labelsize':9,'axes.titlesize':10,'figure.dpi':120,'savefig.dpi':400})
def model(name,group='pooled'):
 q=res[(res.analysis==name)&(res.group==group)];assert len(q)==1,(name,group);return q.iloc[0]
fig,ax=plt.subplots(figsize=(7.0,4.8));rows=[('All 49 items','fi_all','pooled'),('Exclude sleep','fi_no_sleep','pooled'),('Exclude fatigue','fi_no_fatigue','pooled'),('Exclude sleep + fatigue','fi_neither','pooled'),('  2011–2012','fi_neither_cycle','2011-2012'),('  2013–2014','fi_neither_cycle','2013-2014'),('  Additional covariates','fi_neither_fullcase_adjusted','pooled'),('  Exclude mood + care items','fi_no_mood_no_care','pooled'),('  Complete-item cases','fi_neither_complete_cases','pooled'),('  Age ≥65 years','fi_neither_age65','pooled')]
source1=[]
for i,(lab,key,group) in enumerate(rows):
 r=model(key,group);source1.append({'label':lab.strip(),'analysis':key,'group':group,'n':int(r.n),'estimate':r.estimate,'ci_low':r.ci_low,'ci_high':r.ci_high});y=len(rows)-i-1;c='#0072B2' if i<4 else '#555555';ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',ms=4,color=c,ecolor=c,elinewidth=1.2,capsize=2.5);ax.text(.153,y,f"{r.estimate:.3f}  [{r.ci_low:.3f}, {r.ci_high:.3f}]",ha='right',va='center',fontsize=8)
pd.DataFrame(source1).to_csv(P/'figure1_source_data.csv',index=False)
ax.axvline(0,color='#999999',lw=.8);ax.set_yticks(range(len(rows)),[q[0] for q in rows][::-1]);ax.set_ylim(-.8,len(rows)-.15);ax.set_xlim(-.008,.158);ax.set_xticks([0,.025,.05,.075,.10]);ax.set_xlabel('Adjusted difference in mean FI score (0–1)');ax.set_title('Sleep trouble and reconstructed frailty',loc='left',fontweight='bold');ax.text(.153,len(rows)-.2,'Estimate [95% CI]',ha='right',fontsize=8,fontweight='bold');fig.subplots_adjust(left=.33,right=.98,top=.91,bottom=.13)
for ext in ['pdf','svg','png']:fig.savefig(P/f'figure1_fi_associations.{ext}',bbox_inches='tight');plt.close(fig)
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.1,3.35),gridspec_kw={'width_ratios':[1.1,1]});g=d[(d.WTMEC2YR>0)&(d.common_scores==1)];bins=np.arange(0,1.0001,.05)
source2=[]
for group,label,color in [(0,'No','#0072B2'),(1,'Yes','#D55E00')]:
 z=g[g.sleep_told_doctor==group];weights=z.WTMEC2YR.to_numpy()/2;counts,_=np.histogram(z.fi_neither,bins=bins,weights=weights);ax1.stairs(counts/counts.sum(),bins,label=label,color=color,lw=1.6);source2.extend([{'panel':'a','group':label,'bin_lower':bins[j],'bin_upper':bins[j+1],'estimate':counts[j]/counts.sum(),'ci_low':np.nan,'ci_high':np.nan} for j in range(len(counts))])
ax1.set_xlim(0,.85);ax1.set_xlabel('FI excluding sleep and fatigue');ax1.set_ylabel('Weighted share per 0.05 bin');ax1.set_title('a   FI distribution',loc='left',fontweight='bold');ax1.legend(title='Ever told a professional\nabout sleep trouble',frameon=False,fontsize=8,title_fontsize=8)
for i,(key,label,c) in enumerate([('mortality_FI_high_0','≤0.21','#0072B2'),('mortality_FI_high_1','>0.21','#D55E00')]):
 v=des[des.label==key].iloc[0];b=v.estimate*100;lo=v.ci_low*100;hi=v.ci_high*100;source2.append({'panel':'b','group':label,'bin_lower':np.nan,'bin_upper':np.nan,'estimate':v.estimate,'ci_low':v.ci_low,'ci_high':v.ci_high});ax2.errorbar(i,b,yerr=[[b-lo],[hi-b]],fmt='o',color=c,capsize=3,ms=5,lw=1.2);ax2.text(i,b+5,f'{b:.1f}%',ha='center',fontsize=8)
pd.DataFrame(source2).to_csv(P/'figure2_source_data.csv',index=False)
ax2.set_xticks([0,1],['≤0.21','>0.21']);ax2.set_xlim(-.55,1.55);ax2.set_ylim(0,45);ax2.set_ylabel('Weighted mortality status by 2019 (%)');ax2.set_xlabel('Exploratory FI category');ax2.set_title('b   Mortality construct check',loc='left',fontweight='bold');fig.tight_layout(w_pad=2)
for ext in ['pdf','svg','png']:fig.savefig(P/f'figure2_distribution_mortality.{ext}',bbox_inches='tight');plt.close(fig)
print('Wrote vector PDF/SVG and 400 dpi PNG figures')
