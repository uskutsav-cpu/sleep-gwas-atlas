"""Static scientific figures from locked predictions and retained learning curves."""
from common import *
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'runtime_tmp/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    result=json.loads((OUT/'reserved_study_results.json').read_text());rows=list(csv.DictReader((OUT/'reserved_all_observations.csv').open()));figures=OUT/'figures';figures.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220});y=np.array([float(r['measured_absolute_log2_effect']) for r in rows]);cell=np.array([r['cell'] for r in rows]);arms=['sequence_only','AVI_added'];labels=['DNA and cell inputs','AVI added'];colors=['#2166ac','#d95f02'];files=[]
    fig,axes=plt.subplots(2,2,figsize=(9,8),layout='constrained')
    for ci,c in enumerate(['HNP','HEK293']):
        ix=cell==c;limit=max(float(np.log1p(y[ix]).max()),*[float(np.log1p([float(r[a+'_predicted_absolute_log2_effect']) for r in rows if r['cell']==c]).max()) for a in arms])
        for ai,a in enumerate(arms):
            ax=axes[ci,ai];p=np.array([float(r[a+'_predicted_absolute_log2_effect']) for r in rows]);hb=ax.hexbin(np.log1p(y[ix]),np.log1p(p[ix]),gridsize=45,mincnt=1,bins='log',cmap='viridis',extent=(0,limit,0,limit));ax.plot([0,limit],[0,limit],color='gray',ls='--',lw=.8);ax.set(xlim=(0,limit),ylim=(0,limit),xlabel='Measured effect size: log(1 + magnitude)',ylabel='Predicted effect size: log(1 + magnitude)',title=c+' | '+labels[ai]);m=result['arms'][a]['cells'][c];ax.text(.03,.97,f'n = {m["n"]:,}\nMAE = {m["MAE"]:.3f} log2 units\nSpearman = {m["global_Spearman"]:.3f}',transform=ax.transAxes,va='top',bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'});fig.colorbar(hb,ax=ax,label='Measurements per hexagon')
    fig.suptitle('Outside-study effect sizes: all source-eligible measurements',fontsize=12);path=figures/'outside_measured_vs_predicted.png';fig.savefig(path);plt.close(fig);files.append(path)
    fig,axes=plt.subplots(1,3,figsize=(12,3.7),layout='constrained');xx=np.arange(2);width=.34
    for ai,a in enumerate(arms):
        values=[];lower=[];upper=[];coverage=[];accepted=[]
        for c in ['HNP','HEK293']:
            m=result['arms'][a]['cells'][c];values.append(m['MAE']);lo,hi=m['MAE_genomic_window_CI95'];lower.append(m['MAE']-lo);upper.append(hi-m['MAE']);coverage.append(100*m['interval_coverage']);accepted.append(100*m['ranked_fraction'])
        offset=xx+(ai-.5)*width;axes[0].bar(offset,values,width,color=colors[ai],label=labels[ai],yerr=np.maximum([lower,upper],0),capsize=3);axes[1].bar(offset,coverage,width,color=colors[ai]);axes[2].bar(offset,accepted,width,color=colors[ai])
    base=[result['arms']['sequence_only']['cells'][c]['baseline_MAE'] for c in ['HNP','HEK293']];axes[0].scatter(xx,base,color='black',marker='_',s=300,label='Simple cell median',zorder=3);axes[0].set(ylabel='Mean absolute error, log2 units',title='Size-prediction error');axes[0].legend(fontsize=8);axes[1].axhline(90,color='black',ls='--',lw=1);axes[1].set(ylim=(0,105),ylabel='Observed effects inside interval, %',title='Outside interval coverage');axes[2].set(ylim=(0,105),ylabel='Cases eligible for research ranking, %',title='Answers retained by frozen rules')
    for ax in axes:ax.set_xticks(xx,['HNP','HEK293'])
    fig.suptitle('Error, coverage and abstention describe different properties',fontsize=12);path=figures/'outside_error_coverage_and_abstention.png';fig.savefig(path);plt.close(fig);files.append(path)
    fig,ax=plt.subplots(figsize=(9,4),layout='constrained');curve_rows=[]
    recipes=[('cell_effect','Initial magnitude loss','#2166ac'),('raw_MAE','Direct size-error loss','#1b9e77'),('raw_MAE_pairs','Direct loss + measured differences','#d95f02'),('raw_MAE_no_cell_pairs','Generic-cell control','#777777')]
    for recipe,label,color in recipes:
        curves=[]
        for seed in [42,73,101]:
            curve=json.loads((OUT/f'development_{recipe}_{seed}_curve.json').read_text());curves.append(curve)
            for row in curve:curve_rows.append({'recipe':recipe,'seed':seed,**row})
        n=min(len(c) for c in curves);a=np.array([[r['macro_MAE_skill'] for r in c[:n]] for c in curves])*100;epochs=np.arange(1,n+1);ax.plot(epochs,a.mean(0),label=label,color=color);ax.fill_between(epochs,a.min(0),a.max(0),color=color,alpha=.12)
    ax.axhline(0,color='black',lw=.8,ls='--');ax.set(xlabel='Training epoch',ylabel='Size-error improvement over cell median, %',title='Development loci: more training does not ensure better predictions');ax.legend(fontsize=8);path=figures/'development_size_error_learning_curves.png';fig.savefig(path);plt.close(fig);files.append(path)
    columns=sorted(set().union(*(r.keys() for r in curve_rows)));source=figures/'learning_curve_source_data.csv'
    with source.open('w') as f:w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(curve_rows)
    files.append(source);(figures/'figure_notes.txt').write_text('Outside figures use every one of9,602source-eligible measurements from the locked McAfee2023study. Magnitude means absolute source-fitted log2effect; log(1+magnitude) is a display transform only. Full observed ranges are plotted. Shared cells and nearby variants are correlated observations. Error bars resample one-Mb genomic windows, not ancestry-specific LD blocks. The90%line is the nominal calibration target, not a cross-study guarantee. The eligible fraction uses the frozen research-ranking policy and is not a clinical accuracy rate. HNP may be unranked because the development cell group has only one source. Learning-curve shading is the min–max range across three fixed model seeds, using epochs present in all three runs; seeds are not independent biological samples. All retained runs remain in their original curveJSONfiles.\n');files.append(figures/'figure_notes.txt');write_json(OUT/'figure_receipt.json',{'figures_and_sources':{str(p):sha(p) for p in files},'locked_result_sha256':sha(OUT/'reserved_study_results.json'),'joined_source_sha256':sha(OUT/'reserved_all_observations.csv'),'scope':'Static scientific displays of measured performance, no added model fitting.'});print('SCIENTIFIC FIGURES GENERATED',len(files),flush=True)
if __name__=='__main__':main()
