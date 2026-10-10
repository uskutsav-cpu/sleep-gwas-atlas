"""Tiny invented-file tests of the actual scratch disposer, never real spools."""
import copy
import ast
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
REVIEW_PLAN=R/'independent_spool_disposal_plan_v3_1.json'
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
spec=importlib.util.spec_from_file_location('actual_disposal_helper',P/'scripts/86_cleanup_core_bounded_spool_v3.py')
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
          'dependencies_sha256':{str(P/'scripts/86_cleanup_core_bounded_spool_v3.py'):sha(P/'scripts/86_cleanup_core_bounded_spool_v3.py'),str(P/'scripts/extension_replay_common_v3.py'):sha(P/'scripts/extension_replay_common_v3.py')},'archived_input_sha256':{}}
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
    worker={'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':ph,'command':[sys.executable,'-B',str(P/'scripts/85_validate_core_large35_bounded_replay_v3.py'),'--plan',str(path),'--plan-sha',ph,'--trait',trait,'--mode','compare','--out',m['comparison_receipt']], 'returncode':0,'stop_reason':None,'metadata_errors':[],'plan_unchanged':True,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[],'cleanup_error':None,'teardown_verified':True},'output_sha256':{m['comparison_receipt']:sha(m['comparison_receipt'])}}
    save(worker_path,worker)
    return {'folder':folder,'spool':spool,'member':m,'plan':path,'plan_sha':ph,'pieces':pieces,'comparison':comparison,'candidate':candidate,'raw':raw,'old_h':old_h,'old_m':old_m,'worker_path':worker_path,'worker_hash':sha(worker_path),'candidate_hash':sha(spool/'global_harmonization_candidate_receipt.json'),'columns_hash':sha(cp)}
def modify_json(path,callback):
    path=Path(path);d=json.loads(path.read_text());callback(d);path.write_text(json.dumps(d,indent=2)+'\n')
def invoke(f,setup=None,interrupt_after=None,after_inventory=None,mutate_after_first=None,sync_hook=None,expected_output_mutations=(),after_cleanup=None):
    fence();m=f['member'];spool=f['spool'];events=[];deleted=[];fsyncs=[];error=None
    if setup:setup(f)
    before_sources={str(p):sha(p) for p in [f['raw'],f['old_h'],f['old_m'],Path(m['harmonized']),Path(m['munged']),Path(m['harmonization_qc'])]}
    write=h.write_new;unlink=Path.unlink;fsync=common.os.fsync;directory_sync=h.sync_directory;sync_calls=[]
    def durable(fd):
        fsyncs.append({'kind':'directory' if stat.S_ISDIR(os.fstat(fd).st_mode) else 'file'});return fsync(fd)
    def persist(path,data):
        result=write(path,data);events.append('write:'+Path(path).name)
        if str(path)==m['spool_inventory'] and after_inventory:after_inventory(f)
        if str(path)==m['spool_cleanup_receipt'] and after_cleanup:after_cleanup(f)
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

# Mechanical inheritance: apart from declared stricter proof/payload hashing,
# the run body and all other helpers retain the exact predecessor AST.
oldtree=ast.parse((P/'scripts/80_cleanup_core_bounded_spool_v2.py').read_text())
newtree=ast.parse((P/'scripts/86_cleanup_core_bounded_spool_v3.py').read_text())
oldfunc={x.name:x for x in oldtree.body if isinstance(x,ast.FunctionDef)}
newfunc={x.name:x for x in newtree.body if isinstance(x,ast.FunctionDef)}
class Versions(ast.NodeTransformer):
    def visit_Constant(self,node):
        if isinstance(node.value,str):
            node.value=node.value.replace('85_validate_core_large35_bounded_replay_v3.py','79_validate_core_large35_bounded_replay_v2.py').replace('owned_core_ephemeral_spool_inventory_v3','owned_core_ephemeral_spool_inventory_v2').replace('owned_core_ephemeral_spool_cleanup_v3','owned_core_ephemeral_spool_cleanup_v2')
        return node
newfunc=Versions().visit(copy.deepcopy(newtree));newfunc={x.name:x for x in newfunc.body if isinstance(x,ast.FunctionDef)}
for name in oldfunc:
    if name!='run':check('inherited_exact_helper_AST_'+name,ast.dump(oldfunc[name])==ast.dump(newfunc[name]))
def filter_body(body):
    out=[]
    for node in body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in {'retained_hashes','cleanup_hash','intended_inventory_hash'} for t in node.targets):continue
        if isinstance(node,ast.If) and any(isinstance(x,ast.Constant) and isinstance(x.value,str) and x.value in {'NESTED_EXACT_FULL_STREAM_AND_CRC_PROOFS_REQUIRED','PERSISTED_INVENTORY_DIFFERS_FROM_FROZEN_INTENT','PERSISTED_CLEANUP_DIFFERS_FROM_FROZEN_INTENT'} for x in ast.walk(node)):continue
        out.append(node)
    return ast.Module(body=out,type_ignores=[])
check('run_AST_unchanged_except_declared_strict_schema_and_intended_hash_gates',ast.dump(filter_body(oldfunc['run'].body))==ast.dump(filter_body(newfunc['run'].body)))

f=fixture('valid_fixed_intended_inventory_and_copies');r=invoke(f)
inv=json.loads(Path(f['member']['spool_inventory']).read_text());done=json.loads(Path(f['member']['spool_cleanup_receipt']).read_text())
check('ordinary_success_fixed_retained_identities',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'] and done['inventory_sha256']==sha(f['member']['spool_inventory'])==hashlib.sha256((json.dumps(inv,indent=2)+'\n').encode()).hexdigest() and done['copied_candidate_receipt_sha256']==inv['copied_candidate_receipt_sha256']==f['candidate_hash'] and done['copied_column_manifest_sha256']==inv['copied_column_manifest_sha256']==f['columns_hash'] and done['deleted_files']==len(inv['file_inventory']))

rejected('former_inventory_after_write_before_hash_capture',after_inventory=lambda f:modify_json(f['member']['spool_inventory'],lambda d:d.update(file_inventory={},generated_bytes=0)))
def remove_required(f):
    def drop(d):
        for key in ['decompressed_sha256_new','decompressed_sha256_archived','counts_new','counts_archived']:d['munged_comparison'].pop(key)
    corrupted_proof(f,drop)
rejected('former_missing_munged_digest_count_pairs',remove_required)

malformed=[
 ('harmonic_SHA_empty',lambda d:d['harmonized_comparison'].update(decompressed_sha256=['',''])),
 ('harmonic_SHA_nonhex',lambda d:d['harmonized_comparison'].update(decompressed_sha256=['g'*64]*2)),
 ('harmonic_SHA_wrong_type',lambda d:d['harmonized_comparison'].update(decompressed_sha256=[None,None])),
 ('munged_SHA_empty',lambda d:d['munged_comparison'].update(decompressed_sha256_new='',decompressed_sha256_archived='')),
 ('munged_SHA_nonhex',lambda d:d['munged_comparison'].update(decompressed_sha256_new='z'*64,decompressed_sha256_archived='z'*64)),
 ('munged_SHA_short',lambda d:d['munged_comparison'].update(decompressed_sha256_new='a'*63,decompressed_sha256_archived='a'*63)),
 ('munged_SHA_wrong_type',lambda d:d['munged_comparison'].update(decompressed_sha256_new=1,decompressed_sha256_archived=1)),
 ('munged_counts_empty',lambda d:d['munged_comparison'].update(counts_new={},counts_archived={})),
 ('munged_counts_missing_key',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':1},counts_archived={'finite_N_Z':1})),
 ('munged_counts_extra_key',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':1,'missing_or_nonfinite_N_Z':0,'extra':0},counts_archived={'finite_N_Z':1,'missing_or_nonfinite_N_Z':0,'extra':0})),
 ('munged_counts_negative',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':2,'missing_or_nonfinite_N_Z':-1},counts_archived={'finite_N_Z':2,'missing_or_nonfinite_N_Z':-1})),
 ('munged_counts_float',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':1.0,'missing_or_nonfinite_N_Z':0},counts_archived={'finite_N_Z':1.0,'missing_or_nonfinite_N_Z':0})),
 ('munged_counts_bool',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':True,'missing_or_nonfinite_N_Z':0},counts_archived={'finite_N_Z':True,'missing_or_nonfinite_N_Z':0})),
 ('munged_counts_wrong_total',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':2,'missing_or_nonfinite_N_Z':0},counts_archived={'finite_N_Z':2,'missing_or_nonfinite_N_Z':0})),
 ('munged_counts_unequal',lambda d:d['munged_comparison'].update(counts_new={'finite_N_Z':1,'missing_or_nonfinite_N_Z':0},counts_archived={'finite_N_Z':0,'missing_or_nonfinite_N_Z':1}))]
for name,callback in malformed:rejected(name,lambda f,callback=callback:corrupted_proof(f,callback))

f=fixture('cleanup_mutated_immediately_after_write')
r=invoke(f,after_cleanup=lambda f:modify_json(f['member']['spool_cleanup_receipt'],lambda d:d.update(deleted_files=999)))
check('cleanup_intended_hash_rejects_before_first_hash_acceptance','PERSISTED_CLEANUP_DIFFERS_FROM_FROZEN_INTENT' in r['error'] and not r['spool_exists'] and r['success_receipt_exists'])
f=fixture('fixed_inventory_mutated_after_durable_sync')
r=invoke(f,sync_hook=lambda f,count,path:modify_json(f['member']['spool_inventory'],lambda d:d.update(injected=True)) if count==1 else None)
check('fixed_inventory_after_sync_rejects_before_unlink',r['error'] is not None and not r['deleted_files'] and r['spool_exists'])

internal_root=ROOT;ROOT=Path(rp['private_ExFAT_synthetic_root']);ROOT.mkdir(exist_ok=False)
f=fixture('ordinary_actual_ExFAT_transport');sidecars=[str(x.relative_to(f['spool'])) for x in f['spool'].rglob('*') if x.is_file() and x.name.startswith('._')]
r=invoke(f);r['actual_transport_sidecars']=sidecars
check('ordinary_ExFAT_sidecars_durable_sync_and_exact_disposal_pass',r['error'] is None and not r['spool_exists'] and r['success_receipt_exists'] and {'._columns','._post_liftover','._post_ordinary_QC'}.issubset(sidecars) and len(r['directory_sync_calls'])==3)
ROOT=internal_root;fence()
for path,digest in rp['dependencies'].items():check('final_'+Path(path).name,sha(path)==digest)
files={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for base in [ROOT,Path(rp['private_ExFAT_synthetic_root'])] for p in base.rglob('*') if p.is_file() and not p.is_symlink()}
with (R/'independent_spool_disposal_fixture_map_v3_1.json').open('x') as f:json.dump({'synthetic_only':True,'files':files},f,indent=2);f.write('\n')
receipt={'status':'PASS_NARROW_V3_DISPOSAL_CORRECTION_ONLY_NOT_EXECUTION_ADMISSION','checks':CHECKS,'check_count':len(CHECKS),'cases':RESULTS,'case_count':len(RESULTS),'source_large35_plan_sha256':rp['source_large35_plan_sha256'],'plan_sha256':sha(REVIEW_PLAN),'dependencies_before_after':rp['dependencies'],'inherited_v2_33_controls_not_repeated':True,'elapsed_seconds':time.monotonic()-START,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'retained_synthetic_bytes':sum(v['bytes'] for v in files.values()),'production_source_reference_body_reads':0,'production_disposal_mutex_worker_fit_calls':0,'execution_admitted':False}
with (R/'independent_spool_disposal_receipt_v3_1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({k:receipt[k] for k in ['status','check_count','case_count','elapsed_seconds','peak_RSS_bytes','retained_synthetic_bytes']}))
