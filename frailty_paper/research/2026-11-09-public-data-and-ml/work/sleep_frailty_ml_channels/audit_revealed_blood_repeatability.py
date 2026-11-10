"""Repeatability of a previously revealed study, not another model benchmark."""
import itertools,re
from channel_common import *
from scipy.stats import spearmanr,pearsonr
SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'

def read_tsv(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))

def rho(a,b):
    return float(spearmanr(a,b).statistic) if np.ptp(a)>0 and np.ptp(b)>0 else None

def ranges(values):
    values=np.array([v for v in values if v is not None]);return {'minimum':float(values.min()),'median':float(np.median(values)),'maximum':float(values.max()),'splits':len(values)}

def main():
    assert_frozen_previous();assert json.loads((OLD/'validation/reserved_results_independent_audit.json').read_text())['passed'];audit=json.loads((SCREEN/'original_comparator_and_assay_QC_audit.json').read_text());assert audit['reference_and_alternate_source_attributes_reconciled']>0
    cases=list(csv_rows(SCREEN/'candidate_strict_no_published_pretraining_neighborhood.csv'));physical={(r['SNP'],r['comparison'],r['window']):r for r in read_tsv(SCREEN/'Tardaguila2026_sequence_information.tsv')};original={ (r['SNP'],r['comparison'],r['window'],r['Cell_Type']):r for r in read_tsv(SCREEN/'Tardaguila2026_TableS3.tsv')};splits=list(itertools.combinations(range(7),3));assert len(splits)==35;per_cell={};deltas={};export=[];hashes={}
    for cell in ['K562','HL60','CHRF','THP1']:
        count_path=SCREEN/f'{cell}_independent_oligo_count_sums.tsv';counts={r['Oligo']:r for r in read_tsv(count_path)}
        with (SCREEN/f'{cell}_author_sample_metadata.tsv').open() as handle:metadata=list(csv.reader(handle,delimiter='\t'))
        rna={r[0].removesuffix('_cDNA') for r in metadata if r[1]=='RNA'};dna={r[0].removesuffix('_gDNA') for r in metadata if r[1]=='DNA'};assert rna==dna and len(rna)==7;reps=sorted(rna,key=lambda s:(int(re.findall(r'\d+',s)[-1]),s));values=[];targets=[]
        for case in cases:
            key=(case['source_SNP'],case['comparison'],case['window']);design=physical[key];a,b=design['Oligo_ID_REF'],design['Oligo_ID_ALT'];assert a in counts and b in counts;target=float(original[(*key,cell)]['log2Skew']);contrasts=[]
            for rep in reps:
                ra,da,rb,db=[float(counts[allele][rep+suffix]) for allele,suffix in [(a,'_cDNA'),(a,'_gDNA'),(b,'_cDNA'),(b,'_gDNA')]];assert min(ra,da,rb,db)>0
                value=np.log2(rb/db)-np.log2(ra/da);contrasts.append(value);export.append({'source_SNP':key[0],'comparison':key[1],'window':key[2],'cell':cell,'sample_prefix':rep,'original_GLM_signed_coefficient':target,'replicate_log2_ratio_contrast_diagnostic':float(value)})
            values.append(contrasts);targets.append(target)
        delta=np.array(values);target=np.array(targets);deltas[cell]=delta;average=delta.mean(1);signed=[];absolute=[];disagreement=[]
        for split in splits:
            complement=sorted(set(range(7))-set(split));first=delta[:,split].mean(1);second=delta[:,complement].mean(1);signed.append(rho(first,second));absolute.append(rho(abs(first),abs(second)));disagreement.append(float(abs(abs(first)-abs(second)).mean()))
        per_cell[cell]={'constructs':len(cases),'source_rep_names_in_numeric_order':reps,'mean_ratio_vs_original_GLM_signed_Spearman':rho(average,target),'mean_ratio_vs_original_GLM_magnitude_Spearman':rho(abs(average),abs(target)),
                        'mean_ratio_vs_original_GLM_mean_absolute_difference':float(abs(average-target).mean()),'split_half_signed_effect_Spearman':ranges(signed),'split_half_magnitude_Spearman':ranges(absolute),'split_half_magnitude_absolute_disagreement':ranges(disagreement)}
        hashes[str(count_path)]=sha(count_path)
    between={}
    for first,second in itertools.combinations(deltas,2):
        values=[]
        for split in splits:
            complement=sorted(set(range(7))-set(split));a=abs(deltas[first][:,split].mean(1))-abs(deltas[second][:,split].mean(1));b=abs(deltas[first][:,complement].mean(1))-abs(deltas[second][:,complement].mean(1));values.append(rho(a,b))
        between[first+'|'+second]={'split_half_between_cell_magnitude_difference_Spearman':ranges(values),'scope':'Matching numeric positions across independent cell-line replicate lists does not identify biologically paired samples.'}
    from audit_and_lock_library_measurements import write_gzip_csv
    write_gzip_csv(DEST/'revealed_blood_replicate_ratio_diagnostics.csv.gz',export)
    result={'passed':True,'previously_revealed_source_only':True,'original_model_targets_or_results_replaced':False,'new_Ulirsch_outcomes_computed':False,'new_candidate_performance_measured':False,
            'method':'Log2[(alternate RNA/alternate DNA)/(reference RNA/reference DNA)] from original author-filtered counts including their existing per-barcode pseudocount. RNA/DNA library-scale constants cancel within each allele contrast. No new pseudocount or FDR/activity filter. All 35 three-versus-four replicate splits are reported; no favorable split chosen.',
            'per_cell':per_cell,'between_cell_differences':between,'input_count_hashes':hashes,'target_scope':'Simple paired ratio diagnostic, not the original DESeq2 coefficient, a substitute validation target, biological ground truth or a validated ceiling on model accuracy. Split halves still share constructs and barcodes.',
            'code_sha256':sha(Path(__file__))}
    write_json(DEST/'revealed_blood_repeatability_audit.json',result);assert_frozen_previous();print(json.dumps(result),flush=True)

if __name__=='__main__':main()
