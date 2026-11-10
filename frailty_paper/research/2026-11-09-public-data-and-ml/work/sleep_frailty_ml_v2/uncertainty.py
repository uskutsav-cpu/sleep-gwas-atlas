"""Paired chromosome-resampling sensitivity, conditional on fitted selection."""
import json,csv
import numpy as np
from common import *
from sklearn.metrics import roc_auc_score,average_precision_score

def estimate(draws=2000):
    records=json.loads((OUT/'search_results.json').read_text());selected=json.loads((OUT/'model_selection.json').read_text())['selected_model']
    # This comparator is selected on validation, never by external performance.
    dna=[r for r in records if r['features'] in ['kmer','hyena','v1_kmer']]
    comparator=max(dna,key=lambda r:r['validation']['average_precision'])['name']
    names=list(dict.fromkeys([selected,comparator,'avi_direct_rank']));output={}
    datasets=[('reused_benchmark',OUT/'reused_benchmark_predictions.npz',None),('external_assay',WORK/'external_predictions_local.npz',WORK/'external_pairs_local.csv')]
    for label,path,manifest in datasets:
        if not path.exists():continue
        z=np.load(path);y=z['label'];chrom=np.array([int(r['chromosome'].removeprefix('chr')) for r in csv.DictReader(manifest.open())]) if manifest else np.load(OLD/'prepared.npz')['chromosome'][z['row_index']]
        groups=np.unique(chrom);blocks={g:np.flatnonzero(chrom==g) for g in groups};rng=np.random.default_rng(1208);values={n:[] for n in names};skipped=0
        for _ in range(draws):
            ix=np.concatenate([blocks[g] for g in rng.choice(groups,len(groups),replace=True)]);yy=y[ix]
            if len(np.unique(yy))<2:skipped+=1;continue
            for n in names:
                pp=z[n][ix];values[n].append([roc_auc_score(yy,pp),average_precision_score(yy,pp)])
        ranges={n:{metric+'_percentile_range':np.percentile(np.asarray(values[n])[:,j],[2.5,97.5]).tolist() for j,metric in enumerate(['auroc','average_precision'])} for n in names}
        differences={n:{metric+'_percentile_range':np.percentile(np.asarray(values[selected])[:,j]-np.asarray(values[n])[:,j],[2.5,97.5]).tolist() for j,metric in enumerate(['auroc','average_precision'])} for n in names if n!=selected}
        output[label]={'chromosome_clusters':len(groups),'bootstrap_draws':draws,'skipped_draws':skipped,'selected_model':selected,'DNA_only_validation_selected_comparator':comparator,'conditional_model_ranges':ranges,'selected_minus_comparator_ranges':differences,'qualification':'Descriptive dependence sensitivity. Chromosome exchangeability is an assumption; source-assay bias, training, model-selection and foundation-model exposure uncertainty are excluded. Four reused-benchmark groups give particularly unstable ranges.'}
    (OUT/'conditional_uncertainty.json').write_text(json.dumps(output,indent=2)+'\n');return output

if __name__=='__main__':estimate()
