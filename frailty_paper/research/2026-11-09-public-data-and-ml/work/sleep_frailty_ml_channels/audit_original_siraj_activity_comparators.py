"""Check A/B count identities from original tables without assigning GLM signs."""
import re
from channel_common import *
from scipy.stats import spearmanr
SOURCE=DEST/'source_candidates'

def main():
    assert_frozen_previous();receipt=json.loads((SOURCE/'Siraj_original_count_retrieval_receipt.json').read_text())
    for item in receipt['files']:assert sha(item['path'])==item['sha256']
    path=SOURCE/'OL42_20221220_normalized_counts.out.gz';expression={};ambiguous=set();pattern=re.compile(r'(?:^|[;(])([0-9XY]+:[0-9]+:[ACGT]+:[ACGT]+:[RA]:wC)(?=[;):]|$)');n=0
    with gzip.open(path,'rt') as handle:
        reader=csv.reader(handle,delimiter='\t');columns=next(reader);assert len(columns)==20
        for row in reader:
            n+=1;assert len(row)==21;ids=set(pattern.findall(row[0]));values=np.array(row[1:],dtype='float64');assert np.isfinite(values).all() and np.all(values>=0)
            if len(ids)!=1:
                ambiguous.update(ids);continue
            ident=next(iter(ids))
            if ident in expression:ambiguous.add(ident)
            else:expression[ident]=values
    for key in ambiguous:expression.pop(key,None)
    cr=json.loads((V6/'sources/siraj_additional_cell_receipt.json').read_text());files={}
    for f in cr['files']:
        if '_HepG2_' in f['file']:files['HepG2']=(Path(f['file']),f['sha256'],'hepg2')
    neural=BASE/'siraj2026_neural/OL42_SKNSH_emVAR_glm_paired_20221220.out.gz';files['SKNSH']=(neural,sha(neural),'sknsh');results={}
    for cell,(source,expected,prefix) in files.items():
        assert sha(source)==expected;dna=np.array([i for i,k in enumerate(columns) if k.startswith('plasmid_')]);rna=np.array([i for i,k in enumerate(columns) if k.startswith(prefix+'_')]);assert len(dna)==len(rna)==5
        calculated=[];released=[];unavailable=0
        with gzip.open(source,'rt') as handle:
            for r in csv.DictReader(handle,delimiter='\t'):
                if not r['comb'].endswith('_center_fwd_ref'):continue
                rid=r['SNP']+':R:wC';aid=r['SNP']+':A:wC';assert r['ID']==rid
                if rid not in expression or aid not in expression:unavailable+=1;continue
                reference,alternate=expression[rid],expression[aid];calculated.append([reference[dna].mean(),reference[rna].mean(),alternate[dna].mean(),alternate[rna].mean()]);released.append([float(r[k]) for k in ['A_Ctrl_Mean','A_Exp_Mean','B_Ctrl_Mean','B_Exp_Mean']])
        a=np.asarray(calculated);b=np.asarray(released);assert len(a)>100;fields={}
        for i,name in enumerate(['A_control_ref','A_RNA_ref','B_control_alt','B_RNA_alt']):
            difference=abs(a[:,i]-b[:,i]);relative=difference/np.maximum(abs(b[:,i]),1e-8);fields[name]={'same_scale_numerical_matches_at_relative_1e_minus6':int((relative<1e-6).sum()),'median_relative_disagreement':float(np.median(relative)),'maximum_relative_disagreement':float(relative.max()),'rank_correlation':float(spearmanr(a[:,i],b[:,i]).statistic)}
        results[cell]={'original_centered_forward_comparisons_with_unambiguous_single_count_ID_per_allele':len(a),'unavailable_or_multiple_physical_count_IDs':unavailable,'fields':fields,'source_DE_sha256':expected}
    report={'passed_as_count_identity_profile':True,'source_count_rows':n,'unambiguous_center_allele_count_IDs':len(expression),'ambiguous_allele_count_IDs_excluded_from_this_profile':len(ambiguous),'cells':results,
            'no_new_model_fitted_or_target_arrays_changed':True,'original_Log2Skew_sign_convention_established':False,
            'scope':'Original normalized Ref/Alt-labeled count identities compared against A/B means. DESeq2 normalization or filtering may differ between the aggregate count table and per-cell tests. Rank agreement or mean-scale agreement alone does not prove the original paired GLM contrast direction; do not fill unknown signs based on this diagnostic.',
            'code_sha256':sha(Path(__file__))}
    write_json(DEST/'original_siraj_activity_comparator_profile.json',report);assert_frozen_previous();print(json.dumps(report),flush=True)

if __name__=='__main__':main()
