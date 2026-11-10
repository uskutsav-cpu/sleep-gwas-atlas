from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parents[2]/'outputs/sleep_frailty_ai_v1'
r=json.loads((OUT/'test_results.json').read_text());h=json.loads((OUT/'training_history.json').read_text());u=json.loads((OUT/'test_uncertainty.json').read_text())
h.extend({**x,'phase':'small_dna_cnn'} for x in json.loads((OUT/'cnn_training_history.json').read_text()))
models=['kmer_logistic','frozen_hyenadna','finetuned_hyenadna','small_dna_cnn'];names=['Simple DNA\npatterns','Frozen DNA\nencoder','Fine-tuned DNA\nencoder','Local DNA\nconvolutions'];colors=['#4C78A8','#79706E','#E17C05','#59A14F']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(9,4))
for ax,metric,title,ref in zip(axes,['auroc','average_precision'],['Ranking AUC','Average precision'],[.5,r['models']['constant']['prevalence']]):
    vals=[r['models'][m][metric] for m in models]
    ranges=[u['model_ranges'][m][metric+'_percentile_range'] for m in models]
    ax.bar(np.arange(len(models)),vals,color=colors,width=.65)
    for i,(v,q) in enumerate(zip(vals,ranges)):
        ax.vlines(i,q[0],q[1],color='black',lw=1);ax.hlines(q,[i-.07,i-.07],[i+.07,i+.07],color='black',lw=1)
        ax.text(i,max(v,q[1])+.025,f'{v:.3f}',ha='center',fontsize=10)
    ax.axhline(ref,color='#555555',ls='--',lw=1,label='Constant/random ranking')
    ax.set_xticks(range(len(models)),names);ax.tick_params(axis='x',labelsize=8);ax.set_title(title);ax.set_ylim(0,max(.65,max(q[1] for q in ranges)+.10) if metric=='auroc' else max(.12,max(q[1] for q in ranges)+.055));ax.legend(frameon=False,fontsize=8,loc='upper right')
fig.suptitle('Held-out chromosomes 19–22: neuronal reporter calls',fontsize=12)
fig.text(.5,.015,'Whiskers: descriptive chromosome-bootstrap ranges from only four chromosomes; no clinical outcomes.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.05,1,.93])
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'model_comparison.{ext}',dpi=200,bbox_inches='tight')
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(9,3.8))
for phase,color,label in [('frozen_head','#79706E','Head training'),('finetune','#E17C05','Last-block fine-tuning'),('small_dna_cnn','#59A14F','Local DNA CNN')]:
    rows=[x for x in h if x['phase']==phase];xx=[x['epoch'] for x in rows]
    for ax,key in zip(axes,['train_weighted_bce','validation_average_precision']):ax.plot(xx,[x[key] for x in rows],marker='o',ms=3,color=color,label=label)
for ax in axes:ax.set_xlabel('Epoch within phase');ax.legend(frameon=False,fontsize=8)
axes[0].set_ylabel('Training loss (weighted)');axes[1].set_ylabel('Validation average precision')
fig.suptitle('Recorded training curves — test data excluded from selection',fontsize=11);fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'training_curves.{ext}',dpi=200,bbox_inches='tight')
print('Wrote model comparison and training figures')
