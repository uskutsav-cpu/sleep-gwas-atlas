"""Tiny invented-file tests of the actual scratch disposer, never real spools."""
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import resource
import stat
import sys
import time

P=Path(__file__).resolve().parents[1];R=P/'reviews'
REVIEW_PLAN=R/'independent_spool_disposal_plan_v2_1.json'
rp=json.loads(REVIEW_PLAN.read_text());ROOT=Path(rp['private_synthetic_root'])
START=time.monotonic();RESULTS=[];CHECKS=[]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def check(name,condition,detail=None):
    CHECKS.append({'name':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name)
def fence():
    if time.monotonic()-START>rp['bounds']['seconds'] or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>rp['bounds']['RSS_bytes']:raise RuntimeError('REVIEW_BOUND_EXCEEDED')
    if ROOT.exists() and sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file() and not p.is_symlink())>rp['bounds']['namespace_bytes']:raise RuntimeError('REVIEW_NAMESPACE_CAP')
    for path,key in [(P,'internal_floor_bytes'),(ROOT.parent,'SSD_floor_bytes')]:
        s=os.statvfs(path)
        if s.f_bavail*s.f_frsize<rp['bounds'][key]:raise RuntimeError('REVIEW_FREE_FLOOR')
for path,digest in rp['dependencies'].items():check('initial_'+Path(path).name,sha(path)==digest)
fence();ROOT.mkdir(exist_ok=False)
sys.path.insert(0,str(P/'scripts'))
spec=importlib.util.spec_from_file_location('actual_disposal_helper',P/'scripts/80_cleanup_core_bounded_spool_v2.py')
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
import extension_replay_common_v3 as common
def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(data,f,indent=2);f.write('\n')
def gz(path,body):
    path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as f:f.write(body)
def fixture(name):
    folder=ROOT/name;epoch=folder/'private_epoch';trait='synthetic';base=epoch/trait;spool=base/'ephemeral_bounded_spool'
    (spool/'columns').mkdir(parents=True);(spool/'post_liftover').mkdir();(spool/'post_ordinary_QC').mkdir()
    receipts=base/'receipts';receipts.mkdir()
    raw=folder/'synthetic_raw.txt';raw.write_text('invented raw bytes\n')
    old_h=folder/'original_harmonized.gz';new_h=base/'harmonized.gz'
    old_m=folder/'original_munged.gz';new_m=base/'munged.gz'
    hb=b'SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\nrs1\t1\t42\tA\tC\t0.2\t0.1\t0.2\t0.5\t100\n'
    mb=b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t0.500\t100.000\n'
    for path,body in [(old_h,hb),(new_h,hb),(old_m,mb),(new_m,mb)]:gz(path,body)
    qc=base/'QC.tsv';qc.write_text('trait\tsynthetic\nstep\tdropped\tremaining\nall original QC\t0\t1\n')
    m={'trait_id':trait,'ephemeral_spool':str(spool),'harmonized':str(new_h),'munged':str(new_m),'harmonization_qc':str(qc),
       'comparison_receipt':str(receipts/'full_content_QC_comparison.json'),'source_gate_receipt':str(receipts/'source_gate.json'),
       'bounded_receipt':str(receipts/'bounded_harmonization.json'),'column_manifest':str(receipts/'native_column_preparation.json'),
       'spool_inventory':str(receipts/'ephemeral_spool_inventory.json'),'spool_cleanup_receipt':str(receipts/'ephemeral_spool_cleanup.json'),
       'original_design':{'raw':{'sealed_verified_sha256':sha(raw),'resolved_path':str(raw),'bytes':raw.stat().st_size},'harmonizer_input_rows':1,'harmonized_rows':1,
          'harmonized':{'path':str(old_h),'sealed_sha256':sha(old_h)},'munged':{'path':str(old_m),'sealed_sha256':sha(old_m)},'historical_munge_log':{'output_total_rows':1},'historical_QC':{'metadata':{'trait':'synthetic'},'ordered_steps':[{'reason':'all original QC','dropped':0,'remaining':1}]}}}
    m['literal_original_harmonize_command']=[sys.executable,'-B',str(P.parent/'scripts/01_harmonize.py'),'--variant-map','invented_map_not_read']
    plan={'harmonization_python':sys.executable,'private_namespace':str(epoch),'members':[m],'ephemeral_spool_policy':h.POLICY,
          'dependencies_sha256':{str(P/'scripts/80_cleanup_core_bounded_spool_v2.py'):sha(P/'scripts/80_cleanup_core_bounded_spool_v2.py'),str(P/'scripts/extension_replay_common_v3.py'):sha(P/'scripts/extension_replay_common_v3.py')},'archived_input_sha256':{}}
    plan['dependencies_sha256'][str(P.parent/'scripts/01_harmonize.py')]=sha(P.parent/'scripts/01_harmonize.py')
    path=epoch/'synthetic_plan.json';save(path,plan);ph=sha(path)
    pieces=[]
    for relative in ['columns/typed_00_00000000.pkl','columns/tokens_00_00000000.pkl','columns/coordinate_tokens_CHR_00000000.pkl','post_liftover/00000000.pkl','post_ordinary_QC/00000000.pkl']:
        f=spool/relative
        with f.open('xb') as out:pickle.dump({'invented':True,'rows':1},out,protocol=5)
        pieces.append({'path':str(f),'sha256':sha(f),'rows':1})
    (spool/'first_surviving_SNP.sqlite3').write_bytes(b'invented closed database fixture')
    (spool/'required_map_keys.sqlite3').write_bytes(b'invented closed database fixture')
    columns={'expected_rows':1,'source_sha256_before':sha(raw),'source_sha256_after':sha(raw),'typed_pieces':{'SNP':[pieces[0]]},'preinference_token_pieces':{'SNP':[pieces[1]]},'coordinate_token_pieces':{'CHR':[pieces[2]]}}
    cp=spool/'columns/column_preparation_manifest.json';save(cp,columns)
    candidate={'original_code_sha256':sha(P.parent/'scripts/01_harmonize.py'),'status':'CANDIDATE_GLOBAL_HARMONIZATION_PRODUCED_NOT_SCIENTIFICALLY_ADMITTED','source_sha256':sha(raw),'expected_rows':1,'rows_out':1,'output_sha256':sha(new_h),'column_manifest_sha256':sha(cp),'post_liftover_frames':[pieces[3]],'post_ordinary_QC_frames':[pieces[4]],'full_original_output_template_and_QC_comparison_required':True,'original_QC_steps':[['all original QC',0,1]]}
    save(spool/'global_harmonization_candidate_receipt.json',candidate)
    harmonic={'status':'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH','rows':[1,1],'expected_rows':1,'unequal_rows':0,'decompressed_sha256':[hashlib.sha256(hb).hexdigest()]*2,'compressed_sha256':[sha(new_h),sha(old_h)],'both_gzip_CRC_and_EOF_verified':True,'literal_order_fields_missingness_and_numerical_serialization_compared':True,'compressed_byte_identity_required':False}
    munged=common.compare_munged(new_m,old_m,1)
    comparison={'status':'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS','trait':trait,'mode':'compare','plan_sha256':ph,'raw_sha256':sha(raw),'raw_bytes':raw.stat().st_size,'harmonized_comparison':harmonic,'munged_comparison':munged,'QC_comparisons':{'trait':True,'ordered_filter_steps':True},'estimator_calls':0}
    source={'status':'EXACT_CORE_RAW_SOURCE_GATE_PASS','trait':trait,'mode':'source','plan_sha256':ph,'raw_sha256':sha(raw),'raw_bytes':raw.stat().st_size,'estimator_calls':0}
    save(Path(m['source_gate_receipt']),source)
    comparison['output_identity_sha256']={m[k]:sha(m[k]) for k in ['harmonized','munged','harmonization_qc','source_gate_receipt']}
    save(Path(m['comparison_receipt']),comparison)
    worker_path=epoch/'receipts_v4/synthetic__compare.worker.json'
    worker={'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':ph,'command':[sys.executable,'-B',str(P/'scripts/79_validate_core_large35_bounded_replay_v2.py'),'--plan',str(path),'--plan-sha',ph,'--trait',trait,'--mode','compare','--out',m['comparison_receipt']], 'returncode':0,'stop_reason':None,'metadata_errors':[],'plan_unchanged':True,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[],'cleanup_error':None,'teardown_verified':True},'output_sha256':{m['comparison_receipt']:sha(m['comparison_receipt'])}}
    save(worker_path,worker)
    return {'folder':folder,'spool':spool,'member':m,'plan':path,'plan_sha':ph,'pieces':pieces,'comparison':comparison,'candidate':candidate,'raw':raw,'old_h':old_h,'old_m':old_m,'worker_path':worker_path,'worker_hash':sha(worker_path),'candidate_hash':sha(spool/'global_harmonization_candidate_receipt.json'),'columns_hash':sha(cp)}
def modify_json(path,callback):
    path=Path(path);d=json.loads(path.read_text());callback(d);path.write_text(json.dumps(d,indent=2)+'\n')
def invoke(f,setup=None,interrupt_after=None,after_inventory=None,mutate_after_first=None,sync_hook=None,expected_output_mutations=()):
    fence();m=f['member'];spool=f['spool'];events=[];deleted=[];fsyncs=[];error=None
    if setup:setup(f)
    before_sources={str(p):sha(p) for p in [f['raw'],f['old_h'],f['old_m'],Path(m['harmonized']),Path(m['munged']),Path(m['harmonization_qc'])]}
    write=h.write_new;unlink=Path.unlink;fsync=common.os.fsync;directory_sync=h.sync_directory;sync_calls=[]
    def durable(fd):
        fsyncs.append({'kind':'directory' if stat.S_ISDIR(os.fstat(fd).st_mode) else 'file'});return fsync(fd)
    def persist(path,data):
        result=write(path,data);events.append('write:'+Path(path).name)
        if str(path)==m['spool_inventory'] and after_inventory:after_inventory(f)
        return result
    def delete(path,*a,**kw):
        check('only_tiny_exact_spool_unlinked',path.is_relative_to(spool))
        if interrupt_after is not None and len(deleted)==interrupt_after:raise KeyboardInterrupt('synthetic interruption')
        if not deleted:
            check('all_copies_inventory_exist_before_first_unlink',all(Path(m[k]).is_file() for k in ['bounded_receipt','column_manifest','spool_inventory']))
        result=unlink(path,*a,**kw);deleted.append(str(path));events.append('unlink:'+path.relative_to(spool).as_posix())
        if len(deleted)==1 and mutate_after_first:mutate_after_first(f)
        return result
    def sync(path):
        sync_calls.append(str(path));result=directory_sync(path)
        if sync_hook:sync_hook(f,len(sync_calls),path)
        return result
    h.write_new=persist;Path.unlink=delete;common.os.fsync=durable;h.sync_directory=sync
    try:h.run(f['plan'],f['plan_sha'],'synthetic',f['candidate_hash'],f['columns_hash'],f['worker_hash'])
    except BaseException as e:error=type(e).__name__+': '+str(e)
    finally:h.write_new=write;Path.unlink=unlink;common.os.fsync=fsync;h.sync_directory=directory_sync
    check('original_and_current_outputs_preserved_except_declared_fault_injection',all(Path(p).is_file() and (p in expected_output_mutations or sha(p)==v) for p,v in before_sources.items()))
    copies={key:sha(m[key]) for key in ['bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt'] if Path(m[key]).is_file()}
    result={'case':f['folder'].name,'error':error,'deleted_files':deleted,'events':events,'fsync_kinds':fsyncs,'directory_sync_calls':sync_calls,'spool_exists':spool.exists(),'success_receipt_exists':Path(m['spool_cleanup_receipt']).exists(),'retained_copy_hashes':copies}
    RESULTS.append(result);fence();return result

def refresh_worker(f):
    modify_json(f['worker_path'],lambda d:d['output_sha256'].update({f['member']['comparison_receipt']:sha(f['member']['comparison_receipt'])}))
    f['worker_hash']=sha(f['worker_path'])
def corrupted_proof(f,callback):
    modify_json(f['member']['comparison_receipt'],callback);refresh_worker(f)
def rejected(name,setup=None,**hooks):
    f=fixture(name);r=invoke(f,setup,**hooks)
    check(name+'_reject_before_first_unlink',r['error'] is not None and not r['deleted_files'] and r['spool_exists'])
    return f,r

f=fixture('valid_internal_all_paired_sidecars')
for path in [f['spool']/'columns',f['spool']/'post_liftover',f['spool']/'post_ordinary_QC',Path(f['pieces'][0]['path'])]:
    path.with_name('._'+path.name).write_bytes(b'owned invented sidecar')
r=invoke(f);check('valid_internal_pass',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'])
inv=json.loads(Path(f['member']['spool_inventory']).read_text());done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text())
check('valid_retained_byte_hashes_counts_and_namespace',done['deleted_files']==len(inv['file_inventory'])==len(r['deleted_files']) and done['inventory_sha256']==sha(f['member']['spool_inventory']) and done['copied_candidate_receipt_sha256']==inv['copied_candidate_receipt_sha256']==f['candidate_hash'] and done['copied_column_manifest_sha256']==inv['copied_column_manifest_sha256']==f['columns_hash'])
check('real_parent_directory_fsyncs_observed',sum(x['kind']=='directory' for x in r['fsync_kinds'])==3 and len(r['directory_sync_calls'])==3)

for name,setup in [
 ('former_unregistered_typed',lambda f:(f['spool']/'columns/typed_99_99999999.pkl').write_bytes(b'unregistered')),
 ('former_orphan_sidecar',lambda f:(f['spool']/'columns/._typed_99_99999999.pkl').write_bytes(b'orphan')),
 ('former_unregistered_empty_directory',lambda f:(f['spool']/'unregistered_empty').mkdir()),
 ('former_contradictory_scientific_proof',lambda f:corrupted_proof(f,lambda d:(d['harmonized_comparison'].update(status='HARMONIZED_CONTENT_MISMATCH_PRESERVED',unequal_rows=1),d['QC_comparisons'].update(ordered_filter_steps=False)))),
 ('former_stale_munged',lambda f:gz(Path(f['member']['munged']),b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t9.500\t100.000\n')),
 ('former_stale_QC',lambda f:Path(f['member']['harmonization_qc']).write_text('step\tdropped\tremaining\nchanged\t1\t0\n')),
 ('unexpected_plain_file',lambda f:(f['spool']/'unrelated.txt').write_text('preserve')),
 ('unexpected_symlink_file',lambda f:(f['spool']/'columns/typed_98_00000000.pkl').symlink_to(f['raw'])),
 ('unexpected_symlink_directory',lambda f:(f['spool']/'unexpected_link').symlink_to(f['folder'],target_is_directory=True)),
 ('missing_registered_directory',lambda f:(f['spool']/'post_liftover/00000000.pkl').unlink()),
 ('wrong_source_proof',lambda f:modify_json(f['member']['source_gate_receipt'],lambda d:d.update(raw_sha256='0'*64))),
 ('harmonized_output_drift',lambda f:gz(Path(f['member']['harmonized']),b'changed output')),
 ('compare_worker_drift',lambda f:modify_json(f['worker_path'],lambda d:d.update(returncode=1))),
 ('candidate_metadata_drift',lambda f:modify_json(f['spool']/'global_harmonization_candidate_receipt.json',lambda d:d.update(rows_out=2))),
 ('column_metadata_drift',lambda f:modify_json(f['spool']/'columns/column_preparation_manifest.json',lambda d:d.update(expected_rows=2))),
 ('comparison_receipt_drift',lambda f:modify_json(f['member']['comparison_receipt'],lambda d:d.update(trait='other'))),
 ('comparison_failure_addendum',lambda f:Path(str(f['member']['comparison_receipt'])+'.failure.json').write_text('failure')),
 ('CRC_proof_false',lambda f:corrupted_proof(f,lambda d:d['munged_comparison'].update(both_gzip_CRC_and_EOF_verified=False))),
 ('QC_field_missing',lambda f:corrupted_proof(f,lambda d:d['QC_comparisons'].pop('trait')))]:
    rejected(name,setup)
f,r=rejected('former_copied_candidate_mutation_after_inventory',after_inventory=lambda f:modify_json(f['member']['bounded_receipt'],lambda d:d.update(source_sha256='0'*64)))

# The remaining initial-capture gap: modified inventory becomes the frozen hash.
f=fixture('former_inventory_mutation_before_initial_hash_capture')
r=invoke(f,after_inventory=lambda f:modify_json(f['member']['spool_inventory'],lambda d:d.update(file_inventory={},generated_bytes=0)))
inv=json.loads(Path(f['member']['spool_inventory']).read_text());done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text()) if r['success_receipt_exists'] else {}
check('inventory_corruption_before_initial_hash_capture_still_accepted',r['error'] is None and not r['spool_exists'] and done['inventory_sha256']==sha(f['member']['spool_inventory']) and done['deleted_files']>0 and not inv['file_inventory'])

def fail_sync(f,count,path):raise OSError('injected durable parent sync failure')
rejected('durable_parent_sync_failure',sync_hook=fail_sync)
for key in ['bounded_receipt','column_manifest','spool_inventory']:
    def corrupt_at_sync(f,count,path,key=key):
        if count==1:modify_json(f['member'][key],lambda d:d.update(injected_change=True))
    rejected(key+'_change_after_initial_capture_before_unlink',sync_hook=corrupt_at_sync)
for key in ['bounded_receipt','column_manifest','spool_inventory']:
    f=fixture(key+'_change_after_disposal')
    def corrupt_after_disposal(f,count,path,key=key):
        if count==2:modify_json(f['member'][key],lambda d:d.update(injected_change=True))
    r=invoke(f,sync_hook=corrupt_after_disposal)
    check(key+'_post_disposal_change_rejected_with_retained_evidence',r['error'] is not None and not r['spool_exists'] and not r['success_receipt_exists'] and Path(f['member'][key]).exists())
f=fixture('interruption_after_first_unlink');r=invoke(f,interrupt_after=1)
check('interruption_leaves_remaining_scratch_intent_no_success',r['error'].startswith('KeyboardInterrupt') and len(r['deleted_files'])==1 and r['spool_exists'] and not r['success_receipt_exists'] and Path(f['member']['spool_inventory']).exists())
f=fixture('consumed_piece_mutated_during_deletion');r=invoke(f,mutate_after_first=lambda f:Path(f['pieces'][3]['path']).write_bytes(b'changed'))
check('mutated_piece_rejected_during_deletion','GENERATED_SCRATCH_CHANGED' in r['error'] and r['spool_exists'] and not r['success_receipt_exists'])

# A separate schema weakness is isolated from any mutable-hash mismatch.
f=fixture('missing_both_munged_decompressed_digests_and_counts')
def remove_required(f):
    def drop(d):
        for key in ['decompressed_sha256_new','decompressed_sha256_archived','counts_new','counts_archived']:d['munged_comparison'].pop(key)
    corrupted_proof(f,drop)
r=invoke(f,remove_required)
check('missing_required_munged_schema_accepted',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'])

# Real tiny ExFAT transport fixture; all source/data are invented by this checker.
internal_root=ROOT;ROOT=Path(rp['private_ExFAT_synthetic_root']);ROOT.mkdir(exist_ok=False)
f=fixture('valid_actual_ExFAT_transport');before_sidecars=[str(x.relative_to(f['spool'])) for x in f['spool'].rglob('*') if x.is_file() and x.name.startswith('._')]
r=invoke(f);r['actual_transport_sidecars']=before_sidecars
check('actual_ExFAT_sidecars_and_sync_pass',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'] and {'._columns','._post_liftover','._post_ordinary_QC'}.issubset(before_sidecars))
ROOT=internal_root;fence()
for path,digest in rp['dependencies'].items():check('final_'+Path(path).name,sha(path)==digest)
files={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for base in [ROOT,Path(rp['private_ExFAT_synthetic_root'])] for p in base.rglob('*') if p.is_file() and not p.is_symlink()}
with (R/'independent_spool_disposal_fixture_map_v2_1.json').open('x') as f:json.dump({'synthetic_only':True,'files':files},f,indent=2);f.write('\n')
receipt={'status':'BLOCKED_V2_INVENTORY_INITIAL_CAPTURE_AND_NESTED_SCHEMA_GAPS','checks':CHECKS,'check_count':len(CHECKS),'cases':RESULTS,'case_count':len(RESULTS),'source_large35_plan_sha256':rp['source_large35_plan_sha256'],'plan_sha256':sha(REVIEW_PLAN),'dependencies_before_after':rp['dependencies'],'elapsed_seconds':time.monotonic()-START,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'retained_synthetic_bytes':sum(v['bytes'] for v in files.values()),'production_source_reference_body_reads':0,'production_disposal_mutex_worker_fit_calls':0,'execution_admitted':False}
with (R/'independent_spool_disposal_receipt_v2_1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({k:receipt[k] for k in ['status','check_count','case_count','elapsed_seconds','peak_RSS_bytes','retained_synthetic_bytes']}))
