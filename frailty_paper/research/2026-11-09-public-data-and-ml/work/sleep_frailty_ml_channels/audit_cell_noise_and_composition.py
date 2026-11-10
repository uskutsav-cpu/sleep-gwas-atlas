"""Descriptive source-noise diagnosis; no filter, new target, or model selection."""
import itertools
from channel_common import *
from scipy.stats import spearmanr

def profile(y,se):
    keep=np.isfinite(se)&(se>0);a=y[keep];s=se[keep]
    return {'rows':len(y),'positive_original_SE_rows':int(keep.sum()),'mean_absolute_estimate':float(a.mean()),'median_absolute_estimate':float(np.median(a)),
            'mean_original_SE':float(s.mean()),'median_original_SE':float(np.median(s)),
            'mean_absolute_estimate_expected_if_zero_true_effect_and_correct_Gaussian_SE':float(np.sqrt(2/np.pi)*s.mean()),
            'fraction_absolute_estimate_le_SE':float(np.mean(a<=s)),'fraction_absolute_estimate_le_1_96_SE':float(np.mean(a<=1.96*s)),
            'absolute_estimate_SE_Spearman':float(spearmanr(a,s).statistic)}

def main():
    assert_frozen_previous()
    with np.load(OLD/'coverage_training_arrays.npz') as z:d={k:z[k] for k in z.files}
    source=int(np.flatnonzero(d['study_names']=='Siraj')[0]);names=['HepG2','K562','SK-N-SH'];codes={name:list(d['cell_line_names']).index(name) for name in names};rows=np.arange(len(d['effect']));training=(d['split']==0)&(d['study']==source);result={};composition={}
    for library,mask in [('old_libraries',training&(rows<613161)),('added_OL27',training&(rows>=613161))]:
        lookup={line:{int(d['bank_index'][i]):int(i) for i in np.flatnonzero(mask&(d['line']==code))} for line,code in codes.items()};shared=sorted(set.intersection(*[set(x) for x in lookup.values()]));composition[library]={'training_only':True,'shared_three_line_physical_designs':len(shared),'matched_cell_medians':{line:float(np.median(d['effect'][[values[b] for b in shared]])) for line,values in lookup.items()},
                            'matched_mean_cell_differences':{a+'|'+b:float(np.mean(d['effect'][[lookup[a][k] for k in shared]]-d['effect'][[lookup[b][k] for k in shared]])) for a,b in itertools.combinations(names,2)}}
        result[library]={}
        for line,code in codes.items():
            ix=np.flatnonzero(mask&(d['line']==code));assert np.all(d['standard_error_kind'][ix]==2);result[library][line]=profile(d['effect'][ix].astype('float64'),d['measurement_standard_error'][ix].astype('float64'))
    # Already revealed results are used descriptively and never called untouched.
    revealed=list(csv_rows(OLD/'validation/reserved_joined_measurements_and_predictions.csv.gz'));result['revealed_OL33']={}
    for line in names:
        rr=[r for r in revealed if r['source']=='Siraj2026' and r['library']=='OL33' and r['cell_line']==line];result['revealed_OL33'][line]=profile(np.array([float(r['absolute_original_Log2Skew']) for r in rr]),np.array([float(r['original_Skew_SE']) for r in rr]))
    output={'profiles':result,'composition':composition,'original_model_or_partition_changed':False,'new_examples_filtered':False,'SE_used_as_inference_input':False,
            'interpretation':'Cell-wide absolute fitted coefficients can change because true sequence effects, assay conditions and estimation noise all change. The zero-effect expectation assumes unbiased Gaussian coefficient estimates and correctly estimated SE; it is a diagnostic reference, not evidence all effects are zero or a validated denoising target.',
            'scope':'Training-only and explicitly already revealed descriptive audit. Does not replace measured-effect magnitude as the primary objective or justify claims of strong prediction.',
            'input_hashes':{'training':sha(OLD/'coverage_training_arrays.npz'),'revealed':sha(OLD/'validation/reserved_joined_measurements_and_predictions.csv.gz')},'code_sha256':sha(Path(__file__))}
    write_json(DEST/'cell_noise_and_composition_audit.json',output);assert_frozen_previous();print(json.dumps(output),flush=True)

if __name__=='__main__':main()
