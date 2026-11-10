"""Audit public measured neuron activities, without changing the running fits.

This is a source/identity inventory, not a new training dataset authorization.
"""
import collections,csv,gzip,hashlib,json,math,ssl,urllib.request
from pathlib import Path
import numpy as np
import certifi

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_neuronal_activity'
COMMIT='a686f1d552f97ebcc0224057b5a835002f7a3dc8'
FILES=[
 ('cCRE_based_analysis/data/ngn2_element_activity_metadata_202602_bbmap_only_mpralib095_normalized_bcalm_dna_ge_1.tsv.gz',11356492,'aef637e7aba15269a9b30654397d54ba20f70cee'),
 ('cCRE_based_analysis/data/cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_bc10_toptreat_results_scrambled_NP_MK_result_095_dev.tsv.gz',3681474,'d8c2a95634b6b080a5e2d38b43b20e1a2c2905be'),
 ('cCRE_based_analysis/perform_bcalm_elements_cli.R',10767,'ffb8593826c601e508c347b1bd2fe278d11970a0')]


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def rows(path):
    with gzip.open(path,'rt') as handle:yield from csv.DictReader(handle,delimiter='\t' if path.name.endswith('.tsv.gz') else ',')


def digest(sequence):
    reverse=sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]
    return hashlib.sha256(min(sequence,reverse).encode()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True);sources=[];ctx=ssl.create_default_context(cafile=certifi.where())
    for relative,size,expected in FILES:
        path=OUT/Path(relative).name;url=f'https://raw.githubusercontent.com/kircherlab/80K-Analysis/{COMMIT}/{relative}'
        if not path.exists():
            with urllib.request.urlopen(url,context=ctx,timeout=30) as response,path.open('wb') as handle:
                total=0
                for b in iter(lambda:response.read(262144),b''):
                    total+=len(b);assert total<=size;handle.write(b)
        h=hashlib.sha1(f'blob {size}\0'.encode())
        with path.open('rb') as handle:
            for b in iter(lambda:handle.read(1024*1024),b''):h.update(b)
        assert path.stat().st_size==size and h.hexdigest()==expected,path
        sources.append({'path':str(path.relative_to(ROOT)),'url':url,'bytes':size,'Git_blob_SHA1':expected,'sha256':sha(path)})
        print('VERIFIED NEURON ACTIVITY SOURCE',path.name,size,flush=True)
    fitted={}
    for r in rows(OUT/Path(FILES[1][0]).name):
        assert r['name'] not in fitted
        fitted[r['name']]=float(r['logFC'])
    current=ROOT/'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    with np.load(current,allow_pickle=False) as data:
        bank_ids=data['bank_index'];split=data['split'];n=int(bank_ids.max())+1
        masks=np.zeros(n,dtype='uint8');np.bitwise_or.at(masks,bank_ids,(1<<split).astype('uint8'))
        sid=data['study_names'].tolist().index('Salomon');which=data['study']==sid
        known=np.isfinite(data['activity_reference'])&np.isfinite(data['activity_alternate'])
        already={'Salomon_rows_with_existing_measured_activities':int((which&known).sum()),
                 'Salomon_training_rows_with_existing_measured_activities':int((which&known&(split==0)).sum())}
    meta=json.loads((ROOT/'outputs/sleep_frailty_ml_bounded/development/frozen_reader_none_42_started.json').read_text())
    history={};physical_sources=[]
    for source in meta['physical_source_receipt']['sources']:
        path=Path(source['path']);assert sha(path)==source['sha256']
        count=0
        for r in rows(path):
            i=int(r['bank_id'])
            if i>=len(masks) or not masks[i]:continue
            for key in ['reference','alternate']:
                k=digest(r[key]);history[k]=history.get(k,0)|int(masks[i]);count+=1
        physical_sources.append({**source,'allele_sequences_checked':count})
    counts=collections.Counter();lengths=collections.Counter();unique=collections.defaultdict(set)
    activities={};conflicts=set();max_error=0.;output=OUT/'original_neuronal_activity_identity_inventory.csv.gz'
    fields=['name','physical_sequence_sha256_both_orientations','length','source_logFC','history_class','category','chromosome','start','end','strand','source_metadata_sha256']
    with gzip.open(output,'wt') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
        for r in rows(OUT/Path(FILES[0][0]).name):
            counts['metadata_rows']+=1
            try:value=float(r['logFC'])
            except (ValueError,TypeError):counts['missing_activity']+=1;continue
            if not math.isfinite(value):counts['missing_activity']+=1;continue
            assert r['name'] in fitted,r['name'];error=abs(value-fitted[r['name']]);max_error=max(max_error,error);assert error<1e-10
            sequence=r['sequence'].upper()
            if set(sequence)-set('ACGT'):counts['noncanonical_DNA']+=1;continue
            k=digest(sequence);mask=history.get(k,0)
            status='new_physical_sequence' if mask==0 else ('existing_training_only' if mask==1 else 'touches_existing_validation_or_calibration')
            counts[status]+=1;lengths[len(sequence)]+=1;unique[status].add(k)
            if k in activities and abs(activities[k]-value)>1e-6:conflicts.add(k)
            else:activities[k]=value
            writer.writerow({'name':r['name'],'physical_sequence_sha256_both_orientations':k,'length':len(sequence),
               'source_logFC':r['logFC'],'history_class':status,'category':r['category'],'chromosome':r['chr'],
               'start':r['start'],'end':r['end'],'strand':r['strand'],'source_metadata_sha256':sources[0]['sha256']})
    result={'passed_source_and_coefficient_inventory':True,'repository':'kircherlab/80K-Analysis','commit':COMMIT,
        'sources':sources,'existing_target_coverage':already,'counts':dict(counts),'length_counts':dict(lengths),
        'distinct_physical_sequences_by_history':{k:len(v) for k,v in unique.items()},
        'physical_sequences_with_differing_original_construct_coefficients':len(conflicts),
        'max_original_metadata_to_fitted_logFC_difference':max_error,'current_training_arrays_sha256':sha(current),
        'physical_history_sources':physical_sources,'source_target':('Original unshifted BCalm fitted construct log2 RNA/DNA activity coefficient. '
          'The scrambled-control choice determines TREAT statistics; the released logFC is not the difference from the scrambled-control mean. '
          'Not an allelic contrast or phenotype risk.'),
        'targets_inferred_from_AlphaGenome_or_other_model_predictions':False,'source_P_or_activity_threshold_used_by_inventory':False,
        'not_yet_eligible_for_training':'New literal sequences need original genomic-region/component split assignment, validation/calibration neighborhood and whole reserved-study exclusions, duplicate-construct handling, and independent coefficient replay before any new experiment. Current activity and magnitude targets are unchanged.',
        'current_comparison_changed':False,'reserved_Ulirsch_outcome_files_read':False,'new_models_fitted':False,
        'inventory_sha256':sha(output),'code_sha256':sha(Path(__file__))}
    correction=OUT/'activity_target_interpretation_correction.json'
    if correction.exists():
        evidence=json.loads(correction.read_text());assert result['source_target']==evidence['correct_description']
        result['source_target_description_correction']={'receipt':str(correction.relative_to(ROOT)),'sha256':sha(correction),
          'original_receipt_archive':str((OUT/'inventory_receipt_before_activity_target_description_correction.json').relative_to(ROOT))}
    (OUT/'additional_neuronal_activity_inventory_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['counts','distinct_physical_sequences_by_history','max_original_metadata_to_fitted_logFC_difference','physical_sequences_with_differing_original_construct_coefficients']}),flush=True)


if __name__=='__main__':main()
