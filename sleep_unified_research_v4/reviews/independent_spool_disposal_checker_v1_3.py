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
REVIEW_PLAN=R/'independent_spool_disposal_plan_v1_3.json'
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
spec=importlib.util.spec_from_file_location('actual_disposal_helper',P/'scripts/76_cleanup_core_bounded_spool_v1.py')
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
       'comparison_receipt':str(receipts/'compare.json'),'source_gate_receipt':str(receipts/'source.json'),
       'bounded_receipt':str(receipts/'bounded.json'),'column_manifest':str(receipts/'columns.json'),
       'spool_inventory':str(receipts/'inventory.json'),'spool_cleanup_receipt':str(receipts/'cleanup.json'),
       'original_design':{'raw':{'sealed_verified_sha256':sha(raw),'resolved_path':str(raw),'bytes':raw.stat().st_size},'harmonizer_input_rows':1,'harmonized_rows':1,
          'harmonized':{'path':str(old_h),'sealed_sha256':sha(old_h)},'munged':{'path':str(old_m),'sealed_sha256':sha(old_m)},'historical_munge_log':{'output_total_rows':1},'historical_QC':{'metadata':{'trait':'synthetic'},'ordered_steps':[{'reason':'all original QC','dropped':0,'remaining':1}]}}}
    plan={'private_namespace':str(epoch),'members':[m],'ephemeral_spool_policy':h.POLICY,
          'dependencies_sha256':{str(P/'scripts/76_cleanup_core_bounded_spool_v1.py'):sha(P/'scripts/76_cleanup_core_bounded_spool_v1.py'),str(P/'scripts/extension_replay_common_v3.py'):sha(P/'scripts/extension_replay_common_v3.py')},'archived_input_sha256':{}}
    path=folder/'synthetic_plan.json';save(path,plan);ph=sha(path)
    pieces=[]
    for relative in ['columns/typed_00_00000000.pkl','columns/tokens_00_00000000.pkl','columns/coordinate_tokens_CHR_00000000.pkl','post_liftover/00000000.pkl','post_ordinary_QC/00000000.pkl']:
        f=spool/relative
        with f.open('xb') as out:pickle.dump({'invented':True,'rows':1},out,protocol=5)
        pieces.append({'path':str(f),'sha256':sha(f),'rows':1})
    (spool/'first_surviving_SNP.sqlite3').write_bytes(b'invented closed database fixture')
    (spool/'required_map_keys.sqlite3').write_bytes(b'invented closed database fixture')
    columns={'source_sha256_before':sha(raw),'source_sha256_after':sha(raw),'typed_pieces':{'SNP':[pieces[0]]},'preinference_token_pieces':{'SNP':[pieces[1]]},'coordinate_token_pieces':{'CHR':[pieces[2]]}}
    cp=spool/'columns/column_preparation_manifest.json';save(cp,columns)
    candidate={'status':'CANDIDATE_GLOBAL_HARMONIZATION_PRODUCED_NOT_SCIENTIFICALLY_ADMITTED','source_sha256':sha(raw),'expected_rows':1,'rows_out':1,'output_sha256':sha(new_h),'column_manifest_sha256':sha(cp),'post_liftover_frames':[pieces[3]],'post_ordinary_QC_frames':[pieces[4]],'full_original_output_template_and_QC_comparison_required':True,'original_QC_steps':[['all original QC',0,1]]}
    save(spool/'global_harmonization_candidate_receipt.json',candidate)
    harmonic={'status':'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH','rows':[1,1],'expected_rows':1,'unequal_rows':0,'decompressed_sha256':[hashlib.sha256(hb).hexdigest()]*2,'compressed_sha256':[sha(new_h),sha(old_h)],'both_gzip_CRC_and_EOF_verified':True,'literal_order_fields_missingness_and_numerical_serialization_compared':True,'compressed_byte_identity_required':False}
    munged=common.compare_munged(new_m,old_m,1)
    comparison={'status':'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS','trait':trait,'mode':'compare','plan_sha256':ph,'raw_sha256':sha(raw),'raw_bytes':raw.stat().st_size,'harmonized_comparison':harmonic,'munged_comparison':munged,'QC_comparisons':{'trait':True,'ordered_filter_steps':True},'estimator_calls':0}
    source={'status':'EXACT_CORE_RAW_SOURCE_GATE_PASS','trait':trait,'mode':'source','plan_sha256':ph,'raw_sha256':sha(raw),'raw_bytes':raw.stat().st_size,'estimator_calls':0}
    save(Path(m['comparison_receipt']),comparison);save(Path(m['source_gate_receipt']),source)
    return {'folder':folder,'spool':spool,'member':m,'plan':path,'plan_sha':ph,'pieces':pieces,'comparison':comparison,'candidate':candidate,'raw':raw,'old_h':old_h,'old_m':old_m}
def modify_json(path,callback):
    path=Path(path);d=json.loads(path.read_text());callback(d);path.write_text(json.dumps(d,indent=2)+'\n')
def invoke(f,setup=None,interrupt_after=None,after_inventory=None,mutate_after_first=None):
    fence();m=f['member'];spool=f['spool'];events=[];deleted=[];fsyncs=[];error=None
    if setup:setup(f)
    before_sources={str(p):sha(p) for p in [f['raw'],f['old_h'],f['old_m'],Path(m['harmonized']),Path(m['munged']),Path(m['harmonization_qc'])]}
    write=h.write_new;unlink=Path.unlink;fsync=common.os.fsync
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
    h.write_new=persist;Path.unlink=delete;common.os.fsync=durable
    try:h.run(f['plan'],f['plan_sha'],'synthetic')
    except BaseException as e:error=type(e).__name__+': '+str(e)
    finally:h.write_new=write;Path.unlink=unlink;common.os.fsync=fsync
    check('original_and_current_outputs_not_deleted_or_modified',all(Path(p).is_file() and sha(p)==v for p,v in before_sources.items()))
    copies={key:sha(m[key]) for key in ['bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt'] if Path(m[key]).is_file()}
    result={'case':f['folder'].name,'error':error,'deleted_files':deleted,'events':events,'fsync_kinds':fsyncs,'spool_exists':spool.exists(),'success_receipt_exists':Path(m['spool_cleanup_receipt']).exists(),'retained_copy_hashes':copies}
    RESULTS.append(result);fence();return result

f=fixture('full_success_with_exact_partner_sidecar');(f['spool']/'columns/._typed_00_00000000.pkl').write_bytes(b'owned invented AppleDouble')
r=invoke(f);check('full_success',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'])
inv=json.loads(Path(f['member']['spool_inventory']).read_text());done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text())
check('full_success_inventory_counts_and_retained_hashes',done['deleted_files']==len(inv['file_inventory'])==len(r['deleted_files']) and done['inventory_sha256']==sha(f['member']['spool_inventory']) and done['copied_candidate_receipt_sha256']==inv['copied_candidate_receipt_sha256']==sha(f['member']['bounded_receipt']) and done['copied_column_manifest_sha256']==inv['copied_column_manifest_sha256']==sha(f['member']['column_manifest']))
check('actual_helper_has_no_directory_fsync',not any(x['kind']=='directory' for x in r['fsync_kinds']))

for name,setup in [
 ('unregistered_plain_file',lambda f:(f['spool']/'unrelated.txt').write_text('preserve')),
 ('symlink_file',lambda f:(f['spool']/'columns/typed_98_00000000.pkl').symlink_to(f['raw'])),
 ('symlink_directory',lambda f:(f['spool']/'unexpected_dir').symlink_to(f['folder'],target_is_directory=True)),
 ('wrong_source_proof',lambda f:modify_json(f['member']['source_gate_receipt'],lambda d:d.update(raw_sha256='0'*64))),
 ('changed_consumed_piece',lambda f:Path(f['pieces'][0]['path']).write_bytes(b'changed')),
 ('missing_consumed_piece',lambda f:Path(f['pieces'][0]['path']).unlink()),
 ('comparison_top_status_failed',lambda f:modify_json(f['member']['comparison_receipt'],lambda d:d.update(status='FAILED_PRESERVED'))),
 ('candidate_output_hash_wrong',lambda f:modify_json(f['spool']/'global_harmonization_candidate_receipt.json',lambda d:d.update(output_sha256='0'*64)))]:
    f=fixture(name);r=invoke(f,setup)
    check(name+'_reject_preserves_spool',r['error'] is not None and not r['deleted_files'] and r['spool_exists'] and not r['success_receipt_exists'])

# Material counterexamples: each is measured against the actual unchanged helper.
for name,setup in [
 ('unregistered_matching_typed_filename',lambda f:(f['spool']/'columns/typed_99_99999999.pkl').write_bytes(b'unregistered unrelated bytes')),
 ('orphan_matching_AppleDouble',lambda f:(f['spool']/'columns/._typed_99_99999999.pkl').write_bytes(b'orphan unrelated bytes')),
 ('unregistered_empty_directory',lambda f:(f['spool']/'unregistered_empty').mkdir()),
 ('contradictory_nested_scientific_proof',lambda f:modify_json(f['member']['comparison_receipt'],lambda d:(d['harmonized_comparison'].update(status='HARMONIZED_CONTENT_MISMATCH_PRESERVED',unequal_rows=1),d['QC_comparisons'].update(ordered_filter_steps=False)))),
 ('stale_munged_comparison_hash',lambda f:gz(Path(f['member']['munged']),b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t9.500\t100.000\n')),
 ('stale_QC_comparison',lambda f:Path(f['member']['harmonization_qc']).write_text('step\tdropped\tremaining\nchanged QC\t1\t0\n'))]:
    f=fixture(name);r=invoke(f,setup)
    check(name+'_demonstrates_unsafe_acceptance',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'])

f=fixture('copied_candidate_mutated_after_inventory')
r=invoke(f,after_inventory=lambda f:modify_json(f['member']['bounded_receipt'],lambda d:d.update(source_sha256='0'*64)))
inv=json.loads(Path(f['member']['spool_inventory']).read_text());done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text())
check('retained_candidate_mutation_not_detected',r['error'] is None and done['copied_candidate_receipt_sha256']==sha(f['member']['bounded_receipt']) and done['copied_candidate_receipt_sha256']!=inv['copied_candidate_receipt_sha256'])
f=fixture('retained_inventory_mutated_after_write')
r=invoke(f,after_inventory=lambda f:modify_json(f['member']['spool_inventory'],lambda d:d.update(file_inventory={},generated_bytes=0)))
done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text());inv=json.loads(Path(f['member']['spool_inventory']).read_text())
check('retained_inventory_mutation_not_detected',r['error'] is None and done['inventory_sha256']==sha(f['member']['spool_inventory']) and done['deleted_files']>0 and not inv['file_inventory'])

f=fixture('interruption_after_first_unlink');r=invoke(f,interrupt_after=1)
check('interruption_retains_intent_remaining_scratch_no_success',r['error'].startswith('KeyboardInterrupt') and len(r['deleted_files'])==1 and r['spool_exists'] and not r['success_receipt_exists'] and all(Path(f['member'][k]).exists() for k in ['bounded_receipt','column_manifest','spool_inventory']))
f=fixture('piece_changed_during_unlink')
r=invoke(f,mutate_after_first=lambda f:Path(f['pieces'][3]['path']).write_bytes(b'changed during disposal'))
check('mid_disposal_changed_piece_rejected_remaining_proof_preserved','GENERATED_SCRATCH_CHANGED_DURING_DISPOSAL' in r['error'] and r['spool_exists'] and not r['success_receipt_exists'] and Path(f['member']['spool_inventory']).exists())
for path,digest in rp['dependencies'].items():check('final_'+Path(path).name,sha(path)==digest)
files={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in ROOT.rglob('*') if p.is_file() and not p.is_symlink()}
with (R/'independent_spool_disposal_fixture_map_v1_3.json').open('x') as f:json.dump({'synthetic_only':True,'files':files},f,indent=2);f.write('\n')
receipt={'status':'BLOCKED_MATERIAL_DISPOSAL_PROOF_AND_CLOSURE_DEFECTS','check_count':len(CHECKS),'checks':CHECKS,'case_count':len(RESULTS),'cases':RESULTS,'plan_sha256':sha(REVIEW_PLAN),'source_large35_plan_sha256':rp['source_large35_plan_sha256'],'dependencies_before_after':rp['dependencies'],'elapsed_seconds':time.monotonic()-START,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'retained_synthetic_bytes':sum(x['bytes'] for x in files.values()),'production_source_reference_reads':0,'actual_production_spool_disposals':0,'production_mutex_worker_fit_operations':0,'execution_admitted':False}
with (R/'independent_spool_disposal_receipt_v1_3.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({k:receipt[k] for k in ['status','check_count','case_count','elapsed_seconds','peak_RSS_bytes','retained_synthetic_bytes']}))
