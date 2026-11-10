"""Scientific precision–recall curves from the immutable reserved predictions."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from common import *
def main():
    registry=json.loads((OUT/'final_model_registry.json').read_text());report=json.loads((OUT/'reserved_evaluation_report.json').read_text());receipt=json.loads((OUT/'reserved_test_prediction_receipt.json').read_text());assert receipt['prediction_sha256']==sha(OUT/'reserved_test_predictions.npz')==report['prediction_sha256'];exposure=json.loads((OUT/'reserved_exposure_audit.json').read_text());manifest=[]
    for stem in ['salomon','siraj']:
        path=WORK/(stem+'_manifest_local.csv');assert exposure['source_manifest_sha256'][stem]==sha(path);manifest.extend(csv.DictReader(path.open()))
    y=np.array([int(r['label']) for r in manifest]);study=np.array([STUDIES.index(r['study']) for r in manifest]);active=np.array([int(r['active_screen']) for r in manifest]);roles=registry['role_model_names'];detector=roles['overall_assay_call_detector'];variant=roles['within_active_variant_ranker'];records=registry['records'];base_variant=records[variant]['config'].get('base_model',roles['ablation_contrast_no_avi']);call_models=[(detector,'Chosen detector'),(roles['ablation_motif_chromatin_no_avi'],'Matched detector without AVI'),(roles['full_avi_direct'],'AVI score alone'),(roles['pattern_boost_none_depth3'],'Simple DNA patterns')];variant_models=[(variant,'Chosen DNA + AVI model'),(base_variant,'DNA expert alone'),(roles['combination_conditional_avi_only'],'AVI-only conditional model'),(roles['full_avi_direct'],'AVI score alone')]
    with np.load(OUT/'reserved_test_predictions.npz') as z:ix=z['row_index'];pred={k:z[k] for k in z.files if k!='row_index'}
    fig,axes=plt.subplots(2,2,figsize=(12,8.6),sharex=True,sharey=True);colors=['#087e8b','#d26b24','#6552a0','#6f7780'];strict=np.isin(ix,exposure['retained_Siraj_global_row_indices'])
    for s,study_name in enumerate(STUDIES):
        for j,(stratum,suffix,models) in enumerate([('all_calls','__overall',call_models),('measured_active_regions','__active_variant',variant_models)]):
            ax=axes[s,j];mask=(study[ix]==s)&(strict if s==1 else True)&((active[ix]==1) if j else True);truth=y[ix[mask]];used=set()
            for (name,label),color in zip(models,colors):
                if name in used:continue
                used.add(name);p=pred[name+suffix][mask];precision,recall,_=precision_recall_curve(truth,p);ap=report['results'][study_name]['strata'][stratum]['point_metrics'][name]['average_precision'];ax.step(recall,precision,where='post',label=f'{label}  AP={ap:.3f}',color=color,lw=2 if len(used)==1 else 1.3,alpha=.95)
            ax.axhline(truth.mean(),ls='--',lw=1,color='#9fa6ad',label=f'Call rate={truth.mean():.3f}');title=('All tested regions' if not j else 'Already-active regions')+'\n'+('NGN2 · reused study' if s==0 else 'SK-N-SH · strict reserved set');ax.set_title(title,fontsize=11);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.grid(alpha=.15);ax.legend(loc='upper right',fontsize=7.5,frameon=False);ax.spines[['right','top']].set_visible(False)
            if s==1:ax.set_xlabel('Recall: fraction of detected effects recovered')
            if j==0:ax.set_ylabel('Precision: detected-effect fraction')
    fig.suptitle('Locked molecular-assay evaluation',fontsize=16,y=.98);fig.text(.5,.015,'Scores rank laboratory readouts. These are not patient insomnia/frailty predictions or independent genetic-score validation.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.045,1,.95]);png=OUT/'final_precision_recall.png';svg=OUT/'final_precision_recall.svg';fig.savefig(png,dpi=240);fig.savefig(svg);plt.close(fig);(OUT/'final_figure_receipt.json').write_text(json.dumps({'prediction_sha256':report['prediction_sha256'],'registry_sha256':sha(OUT/'final_model_registry.json'),'png_sha256':sha(png),'svg_sha256':sha(svg),'curve_scope':'Selected models plus prespecified matched/relevant comparators. No model reselection from test curves. Other results retained in full metric table.'},indent=2)+'\n');print('Wrote reserved precision-recall figure',flush=True)
if __name__=='__main__':main()
