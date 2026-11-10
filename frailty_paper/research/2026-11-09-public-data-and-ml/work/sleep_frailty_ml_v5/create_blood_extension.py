"""Version the additional K562 arm without mutating the first experiment."""
import shutil,gzip,collections,math,datetime
from common import *
def main():
    target=DATA.parent/'sleep_frailty_ml_v6';target.mkdir(exist_ok=True);out=target/'outputs';out.mkdir(exist_ok=True);work=ROOT/'work/sleep_frailty_ml_v6';work.mkdir(exist_ok=True)
    assert not (OUT/'final_test_evaluation_started.json').exists()
    for p in WORK.iterdir():
        if p.is_file() and p.suffix in ['.py','.R','.txt'] and p.name!='create_blood_extension.py':
            text=p.read_text().replace('sleep_frailty_ml_v5','sleep_frailty_ml_v6').replace('sleep_frailty_ai_v5','sleep_frailty_ai_v6');(work/p.name).write_text(text)
    source_link=target/'sources'
    if not source_link.exists():source_link.symlink_to(SOURCES,target_is_directory=True)
    linked=['sequence_bank.csv.gz','local151_sequences.npz','fixed_feature_bank.npy','kosicki_blind_design.csv','kircher_observations.csv','kircher_sequence_pairs.csv']
    for name in linked:
        p=target/name
        if not p.exists():p.symlink_to(DATA/name)
    for name in ['fixed_features_receipt.json','selected_motifs.json','kosicki_blind_design_audit.json','kircher_preparation_audit.json','full_feature_gradient_check.json']:
        shutil.copyfile(OUT/name,out/name)
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    original=list(csv.DictReader((ROOT/'work/sleep_frailty_ml_v3/siraj_manifest_local.csv').open()));mapping={original[int(r['original_row_index'])]['variant']:r for r in bank if r['sequence_source']=='siraj'}
    rows=list(csv.DictReader((DATA/'measured_observations.csv').open()));fields=list(rows[0]);original_count=len(rows);source=SOURCES/'OL41_42_K562_emVAR_glm_paired_20221220.out.gz';receipt=json.loads((SOURCES/'siraj_pooled_K562_receipt.json').read_text());assert sha(source)==receipt['sha256']
    with gzip.open(source,'rt') as f:raw=[r for r in csv.DictReader(f,delimiter='\t') if r['comb'].endswith('_center_fwd_ref')]
    counts=collections.Counter(r['SNP'] for r in raw);excluded=collections.Counter();added=[]
    for r in raw:
        key='chr'+r['SNP']
        if key not in mapping:excluded['outside_existing_overlap_excluded_sequence_bank']+=1;continue
        if counts[r['SNP']]!=1:excluded['duplicate_variant_in_pooled_readout']+=1;continue
        try:
            aa,bb,ap,bp,skew=(float(r[k]) for k in ['A_log2FC','B_log2FC','A_logPadj_BF','B_logPadj_BF','Log2Skew']);assert np.isfinite([aa,bb,skew]).all() and all(not math.isnan(v) and v>=0 for v in [ap,bp]);active=(abs(aa)>=1 and ap>=2) or (abs(bb)>=1 and bp>=2)
            if active:q=float(r['Skew_logFDR_act']);assert not math.isnan(q) and q>=0;label=int(q>=1)
            else:label=0
        except (ValueError,AssertionError):excluded['unassessable_activity_or_skew']+=1;continue
        b=mapping[key];added.append({'bank_id':b['bank_id'],'study':'Siraj','experiment':'OL41_42_K562','cell':'K562','cell_group':'blood','label':label,'absolute_effect':abs(skew),'activity_mean':.5*(aa+bb),'active_screen':int(active),'pvalue':'','tags':'','parent_locus':'','locus_group':b['locus_group'],'split':b['split']})
    rows.extend(added)
    with (target/'measured_observations.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    audit=json.loads((OUT/'measured_data_audit.json').read_text());summary=collections.defaultdict(lambda:collections.Counter())
    for r in rows:summary[r['cell']][r['split']]+=1;summary[r['cell']][r['split']+'_detected']+=int(r['label'])
    audit.update({'observations':len(rows),'cells':{k:dict(v) for k,v in summary.items()},'additional_K562_exclusions':dict(excluded),'outputs':{n:sha(target/n) for n in ['sequence_bank.csv.gz','local151_sequences.npz','measured_observations.csv']},'extension':'V6 appends measured pooled OL41_42 K562 results using only existing reference-checked, overlap-excluded constructs. Entire original326166row prefix is unchanged; all locus roles unchanged.','original_V5_observations_sha256':sha(DATA/'measured_observations.csv')})
    write_json(out/'measured_data_audit.json',audit)
    write_json(out/'blood_extension_decision.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'The initial exactOL42prefix omitted the pooledOL41_42K562 file. Saturation-only blood rows have no calibration loci in the initial partition; adding the actual matched K562 library enables independent-locus blood calibration rather than unjustified pooled certainty.','original_V5_observations':original_count,'added_blood_observations':len(added),'added_detected_effects':sum(r['label'] for r in added),'total_observations':len(rows),'exclusions':dict(excluded),'source_receipt_sha256':sha(SOURCES/'siraj_pooled_K562_receipt.json'),'old_source_prefix_unchanged':True,'sequence_bank_and_locus_roles_unchanged':True,'reserved_outcomes_opened':False,'old_development_run_preserved':'V5development completes; V5refits intentionally not launched. V6is the expanded-data final arm with the same finite recipe search.','uncertainty_addendum':'Cells lacking>=10loci,>=50rows and>=10ofeachlabel in EACH separate calibration role are explicitly unvalidated/unranked, even if trained. Rare cell types must not inherit an unqualified confidence claim.'})
    link=ROOT/'outputs/sleep_frailty_ai_v6'
    if not link.exists():link.symlink_to(out,target_is_directory=True)
    print('BLOOD EXTENSION CREATED',len(added),'added;',len(rows),'total observations',flush=True)
if __name__=='__main__':main()
