"""Recover original neuronal effect direction without changing current fits."""
from representation_common import *


def main():
    assert_versions_all();folder=DEST/'neural_source_trace';retrieval=json.loads((folder/'retrieval_receipt.json').read_text())
    assert retrieval['passed']
    for r in retrieval['files']:assert sha(r['path'])==r['sha256']
    code=(folder/'perform_bcalm_variants_cli.R').read_text()
    assert 'df_ref$allele <- "ref"' in code and 'df_alt$allele <- "alt"' in code
    assert 'alt = grepl("alt", colnames(mpra))' in code and 'topTable(mpralm_fit, coef = 2, number = Inf)' in code
    with np.load(OLD/'coverage_training_arrays.npz') as z:d={k:z[k] for k in z.files}
    sid=list(d['study_names']).index('Salomon');eligible=(d['study']==sid)&(d['split']!=2)
    obs={int(d['bank_index'][i]):int(i) for i in np.flatnonzero(eligible)};assert len(obs)==int(eligible.sum())
    physical={};pairs=collections.defaultdict(list)
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        bank=int(r['bank_id'])
        if bank in obs:physical[bank]=r;pairs[(r['reference'],r['alternate'])].append(bank)
    assert set(physical)==set(obs)
    result_path=folder/'NGN2_variant_bcalm_80k_bbmap_251128_bcalm_normalized_variant_map_names_dev.tsv'
    original={}
    with result_path.open() as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            assert r['variant_id'] not in original;original[r['variant_id']]=float(r['logFC'])
    mapping={}
    with gzip.open(folder/'2509_variant_map_region_based_controls.tsv.gz','rt') as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            value=(r['REF'],r['ALT']);assert r['ID'] not in mapping or mapping[r['ID']]==value;mapping[r['ID']]=value
    annotation=ROOT/'work/salomon_neuronal_variant_readout.tsv.gz'
    expected={r['path']:r['sha256'] for r in json.loads((V8/'outputs/measured_training_data_audit.json').read_text())['source_files']}
    assert sha(annotation)==expected[str(annotation)]
    matched=collections.defaultdict(dict)
    with gzip.open(annotation,'rt') as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            key=(r['ref_seq'],r['alt_seq'])
            if key not in pairs:continue
            identity=r['ID'];assert identity in original and mapping[identity]==(r['REF'],r['ALT'])
            value=original[identity];assert np.isfinite(value) and abs(value-float(r['NGN2_variant_logFC']))<1e-12
            for bank in pairs[key]:
                index=obs[bank];p=physical[bank];at=int(p['changed_position'])
                assert [i for i,(a,b) in enumerate(zip(*key)) if a!=b]==[at]
                signed=(identity,r['REF'],r['ALT'],value)
                matched[index][signed]={'observation_index':index,'bank_index':bank,'source_variant_id':identity,
                    'source_REF_construct':r['REF'],'source_ALT_construct':r['ALT'],
                    'variant_GRCh38':p['variant_GRCh38'],'signed_ALT_minus_REF_log2_effect':value,'split':int(d['split'][index])}
    rows=[];excluded=collections.Counter();values=np.full(len(d['effect']),np.nan,dtype='float32')
    for index in np.flatnonzero(eligible):
        candidates=matched[int(index)]
        if len(candidates)!=1:excluded['nonunique_or_missing_full_construct_source_identity']+=1;continue
        r=next(iter(candidates.values()));value=r['signed_ALT_minus_REF_log2_effect']
        assert abs(abs(value)-float(d['effect'][index]))<1e-5
        values[index]=value;rows.append(r)
    assert not np.isfinite(values[d['split']==2]).any()
    assert not np.isfinite(d['signed_effect'][np.isfinite(values)]).any()
    np.save(folder/'neural_original_signed_development.npy',values)
    with gzip.open(folder/'neural_original_signed_development.csv.gz','wt',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(folder/'neural_original_direction_recovery_receipt.json',{'passed':True,
       'eligible_existing_development_measurements':int(eligible.sum()),'recovered_original_signed_targets':len(rows),
       'training_rows':int(np.isfinite(values[d['split']==0]).sum()),'validation_rows':int(np.isfinite(values[d['split']==1]).sum()),
       'calibration_rows_added':0,'excluded':dict(excluded),'original_annotated_source_sha256':sha(annotation),
       'original_fitted_results_sha256':sha(result_path),'variant_map_sha256':sha(folder/'2509_variant_map_region_based_controls.tsv.gz'),
       'training_arrays_sha256':sha(OLD/'coverage_training_arrays.npz'),
       'sidecar_sha256':sha(folder/'neural_original_signed_development.npy'),'CSV_sha256':sha(folder/'neural_original_signed_development.csv.gz'),
       'contrast':'Original fitted log2 RNA/DNA coefficient for ALT indicator, relative to REF.',
       'original_fitted_coefficients_match_existing_annotation_without_sign_transformation':True,
       'source_REF_ALT_constructs_and_full_DNA_match_existing_model_inputs':True,
       'selected_by_significance_activity_or_effect_size':False,'original_counts_or_fitting_reconstructed':False,
       'exact_historical_BCalm_runtime_bound':False,'current_diagnostic_models_used_new_neural_sidecar':False,
       'all_previous_targets_partitions_models_and_results_preserved':True,'new_reserved_test_outcomes_read':False})
    assert_versions_all();print('ORIGINAL NEURONAL SIGNED TARGETS RECOVERED',len(rows),dict(excluded),flush=True)


if __name__=='__main__':main()
