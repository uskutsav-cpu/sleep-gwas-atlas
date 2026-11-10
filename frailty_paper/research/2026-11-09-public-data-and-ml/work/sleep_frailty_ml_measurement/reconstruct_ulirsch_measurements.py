"""Sealed count-based reanalysis; never a claim of recovering original coefficients."""
import re,time
from measurement_common import *
from scipy.special import polygamma

def main():
    assert_prior_versions();started=time.monotonic();out=DEST/'Ulirsch_sealed_count_measurements.csv.gz';assert not out.exists()
    source_audit=json.loads((SOURCE/'Ulirsch_complete_history_and_pretraining_guard.json').read_text());assert source_audit['passed'];cases=list(csv_rows(SOURCE/'Ulirsch_candidate_strict_no_history_or_published_pretraining.csv'));assert len(cases)==source_audit['strict_candidate_designs']
    plan={'source':'Ulirsch et al. 2016, original GEO GSE87711 raw processed counts; new analysis, not original published fitted effects.',
          'physical_and_exposure_source_sha256':sha(SOURCE/'Ulirsch_candidate_strict_no_history_or_published_pretraining.csv'),
          'primary_condition':'Basal K562 only: six CTRL RNA preparations and two DNA preparations. GATA1 overexpression is excluded from this primary reconstruction, not relabeled ordinary K562.',
          'aggregation':'Sum the 14 source count rows within each exact source allele construct, preserving every RNA and DNA sample. Match chr,pos,ref,alt,type,bot,top,clean and actual design ID; never join on the ambiguous byallele field alone.',
          'primary_log2_contrast':'Mean of six log2((ALT_RNA+0.5)/(REF_RNA+0.5)) minus mean of two log2((ALT_DNA+0.5)/(REF_DNA+0.5)). Map source Ref/Mut to the full-verified genomic REF/ALT; magnitude is absolute contrast.',
          'pseudocount':.5,'count_availability':'Both allele constructs must have all 14 original rows, and each allele must have >=10 total DNA counts in each of two DNA samples. At least one RNA count across both alleles/six RNA samples. No significance, activity, measured effect, estimated SE or model-prediction cutoff.',
          'zero_RNA_allele':'Retain a single completely unobserved RNA allele as a regularized, flagged estimate, not a precise finite biological effect. Both RNA alleles absent means contrast unavailable, not zero effect.',
          'empirical_variance':'Sample variance of six RNA log2 allelic ratios divided by six plus sample variance of two DNA log2 allelic ratios divided by two; independent-preparation working assumption.',
          'count_variance_floor':'Sum of trigamma(count+0.5) for both alleles over all RNA samples / (6^2*ln(2)^2), plus corresponding DNA sum / (2^2*ln(2)^2). A Jeffreys count-model log-ratio reference, not a source GLM SE or guaranteed coverage.',
          'reported_SE':'sqrt(max(empirical_variance,count_variance_floor)); all components exported separately. No p value, discovery hit label or original GLM sign assigned.',
          'interpretation':'Observed, regularized RNA/DNA count contrast with transparent uncertainty diagnostics. It is not denoised biological ground truth and cannot validate clinical scores or causal genes.',
          'outcome_policy':'Save coefficient file sealed from model-selection scripts; disclose only availability and reconstruction agreement. Any future test needs frozen candidate, confidence and predictions before joining these values.',
          'new_model_predictions_or_performance_used':False,'code_sha256':sha(Path(__file__))}
    plan_path=DEST/'Ulirsch_count_reanalysis_protocol.json'
    if plan_path.exists():assert json.loads(plan_path.read_text())==plan
    else:write_json(plan_path,plan)
    original=json.loads((SOURCE/'Ulirsch_original_retrieval_receipt.json').read_text())
    for item in original['files']:assert sha(item['path'])==item['sha256']
    design={};pattern=re.compile(r'^([^:]+):(\d+):([ACGT]+):([ACGT]+)_(Ref|Mut)_Bot_([^_]+)_Top_([^_]+)(_mod)?$')
    with gzip.open(SOURCE/'GSE87711_design.txt.gz','rt') as handle:
        for identifier,sequence in csv.reader(handle,delimiter='\t'):
            match=pattern.fullmatch(identifier);assert match;c,p,ra,aa,kind,bot,top,mod=match.groups();key=(c,p,ra,aa,kind,bot,top,'mod' if mod else 'var');assert key not in design;design[key]=identifier
    sample_names=['K562_minP_DNA1','K562_minP_DNA2']+[f'K562_CTRL_minP_RNA{i}' for i in range(1,7)];sums={};nrows=collections.Counter()
    with gzip.open(SOURCE/'GSE87711_RBC_MPRA_minP_raw.txt.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=tuple(row[k] for k in ['chr','pos','ref','alt','type','bot','top','clean']);identifier=design[key];values=np.array([int(row[k]) for k in sample_names],dtype='int64');assert np.all(values>=0);sums.setdefault(identifier,np.zeros(8,dtype='int64'));sums[identifier]+=values;nrows[identifier]+=1
    output=[];availability=collections.Counter()
    for i,case in enumerate(cases):
        ref_id,alt_id=case['source_reference_construct'],case['source_alternate_construct'];assert nrows[ref_id]==nrows[alt_id]==14;reference,alternate=sums[ref_id],sums[alt_id]
        if case['GRCh38_reference_is_source_Ref']!='True':reference,alternate=alternate,reference
        dna_ok=bool(np.all(reference[:2]>=10) and np.all(alternate[:2]>=10));rna_ok=int(reference[2:].sum()+alternate[2:].sum())>0;eligible=dna_ok and rna_ok;availability['DNA_depth_eligible']+=dna_ok;availability['any_RNA_available']+=rna_ok;availability['eligible_point_contrasts']+=eligible
        beta=se=empirical=countvar=float('nan');rna_zero=bool(reference[2:].sum()==0 or alternate[2:].sum()==0)
        if eligible:
            ratios=np.log2((alternate+.5)/(reference+.5));beta=float(ratios[2:].mean()-ratios[:2].mean());empirical=float(ratios[2:].var(ddof=1)/6+ratios[:2].var(ddof=1)/2)
            variance=polygamma(1,alternate+.5)+polygamma(1,reference+.5);countvar=float((variance[2:].sum()/36+variance[:2].sum()/4)/np.log(2)**2);se=float(np.sqrt(max(empirical,countvar)));assert np.isfinite([beta,se]).all() and se>0
        availability['one_RNA_allele_completely_unobserved_retained']+=int(eligible and rna_zero)
        output.append({'design_index':i,'variant_GRCh38':case['variant_GRCh38'],'original_GRCh37_variant':case['original_GRCh37_variant'],'window_bot':case['window_bot'],'window_top':case['window_top'],'component':case['joint_component_after_all_source_DNA'],
                       'source_reference_construct':ref_id,'source_alternate_construct':alt_id,'GRCh38_reference_is_source_Ref':case['GRCh38_reference_is_source_Ref'],'eligible_point_contrast':eligible,'minimum_DNA_count_across_two_alleles_and_samples':int(min(reference[:2].min(),alternate[:2].min())),
                       'at_least_one_RNA_count':rna_ok,'one_RNA_allele_completely_unobserved':rna_zero,'signed_genomic_ALT_minus_REF_log2_ratio':beta,'absolute_log2_ratio':abs(beta),'working_SE':se,'empirical_variance':empirical,'count_variance_floor':countvar,
                       **{f'genomic_REF_{name}':int(value) for name,value in zip(sample_names,reference)},**{f'genomic_ALT_{name}':int(value) for name,value in zip(sample_names,alternate)}})
    with gzip.open(out,'wt',newline='') as handle:writer=csv.DictWriter(handle,fieldnames=list(output[0]));writer.writeheader();writer.writerows(output)
    write_json(DEST/'Ulirsch_sealed_reconstruction_receipt.json',{'passed':True,'candidate_designs':len(cases),'availability':dict(availability),'outcome_file_sha256':sha(out),'protocol_sha256':sha(plan_path),
              'original_published_coefficient_reproduced':False,'new_model_performance_evaluated':False,'new_outcomes_used_for_model_selection':False,'outcome_values_or_distributions_printed':False,'independent_R_check_pending':True,'elapsed_seconds':time.monotonic()-started})
    assert_prior_versions();print('SEALED COUNT RECONSTRUCTION COMPLETE',json.dumps({'designs':len(cases),'availability':dict(availability),'outcome_values_printed':False}),flush=True)

if __name__=='__main__':main()
