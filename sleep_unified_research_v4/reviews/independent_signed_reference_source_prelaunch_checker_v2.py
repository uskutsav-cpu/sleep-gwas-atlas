"""Bounded, invented v2 compact/controller delta review; no real bodies."""
import ast
from collections import Counter
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import struct
import sys
import time
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
PLAN=P/'manifests/signed_reference_source_execution_plan_v2.json'
PLAN_SHA='49b8fed7a8af638362dbbededbb9a01e6877f9b727655a3b7327fdfd07c0ccee'
FIX=P.parents[1]/'independent_signed_reference_source_fixtures_v2'
LIMIT=json.loads((R/'independent_signed_reference_source_prelaunch_resource_plan_v2.json').read_text())
sys.path.insert(0,str(S))
import signed_reference_compact_v2 as compact
import signed_reference_source_common_v2 as common
import signed_reference_source_controller_v2 as controller
import signed_reference_source_worker_v2 as worker
import terminal_commit_common_v2 as terminal

START=time.monotonic();CASES=[];ASSERTIONS=0


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()


def ck(test,message):
    global ASSERTIONS
    ASSERTIONS+=1
    if not test:raise AssertionError(message)


def guard():
    ck(time.monotonic()-START<LIMIT['deadline_seconds'],'review deadline')
    ck(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<LIMIT['RSS_stop_bytes'],'review RSS')
    ck(sum(p.stat().st_size for p in FIX.rglob('*') if p.is_file())<LIMIT['fixture_cap_bytes'],'review storage')


def case(name,fn):
    before=ASSERTIONS;detail=fn();guard()
    CASES.append(dict(case=name,status='PASS',assertions=ASSERTIONS-before,detail=detail))


def put(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
    return sha(p)


def canonical(row):
    return (json.dumps(row,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


def funcs(path):
    return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(path).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}


def unchanged():
    ck(sha(PLAN)==PLAN_SHA,'actual plan')
    p=json.loads(PLAN.read_text());fresh={}
    for path,h in p['dependency_sha256'].items():
        if path not in p['preparation_deferred_dependency_sha256']:
            ck(sha(path)==h,'plan nonbody dependency '+path);fresh[path]=h
    old=json.loads((P/'manifests/signed_reference_source_execution_plan_v1.json').read_text())
    for field in ['candidate','archive_members','reference_only_QC','independent_controls','finite_reference_truth_contract','preserved_full_protocol','reference_body_sha256','expected_master_reference_rows','source_python','decoderB_python','runtime_B_profile_sha256']:
        ck(p[field]==old[field],'unchanged scientific plan '+field)
    for base,names in [('signed_reference_source_controller',['cleanup','runtime_profile_gate','monitor_worker','monitor_helper','group_exists']),('signed_reference_source_worker',['line','consume_jsonl','acquire','configure_db','source_format','read_dose']),('signed_reference_ldsc_decoderB',['export','main'])]:
        a=funcs(S/(base+'_v1.py'));b=funcs(S/(base+'_v2.py'))
        for name in names:ck(a[name]==b[name],'unchanged AST '+base+':'+name)
    ck(sha(S/'signed_reference_decode_v1.py')==sha(S/'signed_reference_decode_v2.py'),'decoder/QC/operators byte-identical')
    oldseal=R/'independent_signed_reference_source_prelaunch_review_v1_seal.json'
    ck(sha(oldseal)=='4b024a5f1583b5d7db08f1b20b0cad6f211946a83dde9788d93af6f501b78e8b','prior rejection seal')
    for path,h in json.loads(oldseal.read_text())['file_sha256'].items():ck(sha(path)==h,'prior rejected artifact unchanged')
    ck(sha(R/'independent_decoderB_runtime_identity_review_v1_seal.json')=='7d7f388d4380632682f3dde920b32f93e2b7e454d90b2a72dac1afde3c7a0a51','completed B identity review reused')
    return dict(current_nonbody_dependency_count=len(fresh),current_nonbody_dependency_sha256=fresh,deferred_body_binary_count=len(p['preparation_deferred_dependency_sha256']),old_suites_and_runtime_census_repeated=0)


def row(rank,chrom,number,reason=None):
    v=dict(SNP=f'00_INVENTED_é_χ_{chrom}_{number}',CHR=chrom,BP=100001+number,master_rank0=rank,master_rank1=rank+1)
    if reason:v['exclusion']=reason
    else:
        v.update(target_A1='A',target_A2='C',source_A1='A' if number%2==0 else 'C',source_A2='C' if number%2==0 else 'A',sign=1 if number%2==0 else -1,source_row0=number,BED_byte_offset=3+126*number,source_member=f'1000G.EUR.QC.{chrom}.bed')
    return v


def q(obs):
    n=len(obs);total=sum(obs);sq=sum(d*d for d in obs);num=n*sq-total*total;miss=503-n
    reason='ALL_MISSING' if not n else 'FEWER_THAN_TWO_OBSERVED' if n<2 else 'EXACT_ZERO_VARIANCE' if num<=0 else 'MISSING_GT25_OF503' if miss>25 else None
    f=total/(2*n) if n else None
    return dict(observed_count=n,missing_count=miss,sum_dose=total,sum_dose2=sq,integer_variance_numerator=num,A1_frequency=f,MAF=min(f,1-f) if f is not None else None,exclusion=reason)


MASTER_ROWS=[];QC_ROWS=[];BASE=FIX/'complete_invented_compact';SB=None;QB=None


def reconstruction():
    global SB,QB
    BASE.mkdir();rank=0;idx=0
    normal=[i%3 for i in range(503)]
    for chrom in range(1,23):
        for i in range(192):
            r=row(rank,chrom,i);MASTER_ROWS.append(r);rank+=1
            item={**r,**q(normal),'J_index0':idx};idx+=1;QC_ROWS.append(item)
    specials=[[],[1],[0]*503,normal[:477],normal[:478],[1]+[0]*502]
    for i,obs in enumerate(specials):
        r=row(rank,22,192+i);MASTER_ROWS.append(r);rank+=1;item={**r,**q(obs)}
        if item['exclusion'] is None:item['J_index0']=idx;idx+=1
        QC_ROWS.append(item)
    for i,reason in enumerate(compact.STATIC_REASONS[1:]):MASTER_ROWS.append(row(rank,22,198+i,reason));rank+=1
    sw=compact.StaticWriter(BASE,len(MASTER_ROWS))
    for r in MASTER_ROWS:sw.append(r)
    SB=sw.close();gw=compact.GenotypeWriter(BASE,SB)
    for r in QC_ROWS:gw.append(r)
    QB=gw.close()
    expected=[([r for r in MASTER_ROWS if 'exclusion' not in r],list(compact.static_rows(BASE,SB)),'static_selection'),([r for r in MASTER_ROWS if 'exclusion' in r],list(compact.static_rows(BASE,SB,selected=False)),'static_exclusions'),([r for r in QC_ROWS if r['exclusion'] is None],list(compact.genotype_rows(BASE,SB,QB)),'J_ordered'),([r for r in QC_ROWS if r['exclusion'] is not None],list(compact.genotype_rows(BASE,SB,QB,eligible=False)),'genotype_exclusions')]
    hashes={}
    for want,got,name in expected:
        ck(want==got,'independent complete dictionary reconstruction '+name)
        hashes[name]=hashlib.sha256(b''.join(canonical(r) for r in want)).hexdigest()
        ck(hashes[name]==(SB if name.startswith('static') else QB)['logical_sha256'][name],'independent canonical v1 bytes '+name)
        for a,b in zip(want,got):
            if name in ['J_ordered','genotype_exclusions']:
                for k in ['A1_frequency','MAF']:
                    ck((a[k] is None and b[k] is None) or (a[k] is not None and b[k] is not None and struct.pack('<d',a[k])==struct.pack('<d',b[k])),'binary64 AF/MAF bit equality')
    ck(SB['master_rows']==4238 and SB['static_rows']==QB['QC_rows']==4230 and QB['J_rows']==4226,'complete counts/order')
    # Independent fixed-format parser, not candidate read routines.
    ownM=struct.Struct('<BQQIQb4sB');ownQ=struct.Struct('<IHHHIB');ownH=struct.Struct('<8sQ')
    mb=(BASE/'master_static.bin').read_bytes();pool=(BASE/'master_SNP_utf8.bin').read_bytes();gb=(BASE/'genotype_QC.bin').read_bytes()
    ck(ownH.unpack(mb[:16])==(b'SRMSTR02',4238) and ownH.unpack(gb[:16])==(b'SRGENQ02',4230),'independent headers')
    for i,r in enumerate(MASTER_ROWS):
        c,bp,off,n,src,sign,alleles,code=ownM.unpack_from(mb,16+35*i)
        ck((c,bp,pool[off:off+n].decode('utf8'))==(r['CHR'],r['BP'],r['SNP']),'independent physical master')
        ck(code==compact.STATIC_REASONS.index(r.get('exclusion')),'independent reason enum')
    for i,r in enumerate(QC_ROWS):
        values=ownQ.unpack_from(gb,16+15*i)
        ck(values==(r['master_rank0'],r['observed_count'],r['sum_dose'],r['sum_dose2'],r['integer_variance_numerator'],compact.QC_REASONS.index(r['exclusion'])),'independent physical exact QC')
    return dict(complete_master_rows=4238,complete_genotype_rows=4230,J_rows=4226,logical_sha256=hashes,binary64_bits_exact=True,actual_source_rows=0)


def bounds_and_failures():
    fixed=1290028*(35+15)+32
    total=fixed+(64+128+64+32+32+16)*(1<<20)
    ck(fixed==64501432 and total==416822968 and 512*(1<<20)-total==120047944,'independent structural arithmetic')
    p=FIX/'pool_fail';p.mkdir();sw=compact.StaticWriter(p,1);sw.pool.cap=2
    try:
        try:sw.append(row(0,1,0))
        except RuntimeError as e:ck('CAP' in str(e),'pool cap whole-stage rejection')
        else:raise AssertionError('pool cap refusal required')
    finally:sw.abort()
    ck((p/'master_static.bin').stat().st_size==16 and not (p/'genotype_QC.bin').exists(),'before QC, no truncated master acceptance')
    # Explicit malformed rank/count values must fail, rather than renumber/drop.
    p=FIX/'rank_fail';p.mkdir();sw=compact.StaticWriter(p,1)
    try:
        r=row(1,1,0)
        try:sw.append(r)
        except RuntimeError as e:ck('RANK' in str(e),'original rank preserved')
        else:raise AssertionError('rank refusal required')
    finally:sw.abort()
    return dict(fixed_bytes=fixed,partition_bound=total,strict_QC_cap=512*(1<<20),structural_margin=120047944,source_dependent_pool_DB_fit_not_observed=True)


def curl_command():
    seen=[]
    class Stop(Exception):pass
    def capture(command,**kw):seen.append(command);raise Stop()
    with patch.object(worker.subprocess,'run',capture):
        try:worker.curl_file({'curl':'/usr/bin/curl','candidate':{'bytes':288277344}},'https://example.invalid/archive',FIX/'not_created')
        except Stop:pass
    ck(seen[0][1]=='-q','first-position q')
    ck(seen[0][seen[0].index('--proto')+1]=='=https','unchanged HTTPS')
    ck(not any(v in seen[0] for v in ['--retry','--location','--insecure','--continue-at']),'no altered science/source transport scope')
    return dict(captured_argv=seen[0],network_operations=0,existing_v1_actual_curl_config_proof_reused=True)


@contextmanager
def no_lock(*args,**kwargs):yield 91


def full_controller(name,fault=None):
    outer=FIX/name;outer.mkdir();root=outer/'canonical_signed_reference_source_v2'
    lock=outer/'INVENTED_NOT_FLOCK';lock.write_bytes(b'INVENTED\n')
    plan=dict(output_root=str(root),expected_master_reference_rows=len(MASTER_ROWS),source_python='INVENTED_PYTHON',source_worker=str(S/'signed_reference_source_worker_v2.py'),candidate=dict(bytes=14,upstream_md5='md5:INVENTED',file_name='INVENTED'),resource_contract=dict(shared_heavy_mutex=str(lock)),raw_transfer_family_mutex=str(lock),compact_representation=compact.SCHEMA)
    pp=outer/'plan.json';ap=outer/'admission.json';ps=put(pp,plan);ads=put(ap,{'INVENTED_NO_ADMISSION':True});initial=[]
    counts=dict(Counter(str(r['CHR']) for r in QC_ROWS if r['exclusion'] is None))
    def fake_monitor(p,command,fds,started,ownership,gates):
        initial.extend(command)
        (root/'extracted/members').mkdir();(root/'qc/controls').mkdir()
        for target in controller.expected_worker_outputs(root):Path(target).write_bytes(b'INVENTED_ONLY\n')
        for f in ['master_static.bin','master_SNP_utf8.bin','genotype_QC.bin']:shutil.copyfile(BASE/f,root/'qc'/f)
        def phase(name,value): (root/'qc'/name).write_text(json.dumps(value,indent=2)+'\n')
        phase('compact_schema.json',compact.SCHEMA)
        ai=dict(bytes=14,md5='INVENTED',sha256=sha(root/'archive/1000G_Phase3_plinkfiles.tgz'))
        phase('G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',candidate=plan['candidate'],archive_identity=ai))
        phase('G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,chromosome_variant_counts=counts,source_index_sha256=sha(root/'extracted/source_index.sqlite')))
        phase('G4_static_seal.json',dict(status='FROZEN_BEFORE_GENOTYPE_QC',master_reference_rows=len(MASTER_ROWS),static_rows=SB['static_rows'],compact_static=SB,compact_schema_file_sha256=sha(root/'qc/compact_schema.json'),exclusion_counts=dict(Counter(r['exclusion'] for r in MASTER_ROWS if 'exclusion' in r)),static_ordered_sha256=SB['logical_sha256']['static_selection'],static_exclusions_sha256=SB['logical_sha256']['static_exclusions']))
        phase('G4_J_seal.json',dict(status='REFERENCE_ONLY_J_FROZEN',compact_genotype=QB,J_count=QB['J_rows'],J_ordered_sha256=QB['logical_sha256']['J_ordered'],genotype_exclusions_sha256=QB['logical_sha256']['genotype_exclusions'],genotype_exclusion_counts=dict(Counter(r['exclusion'] for r in QC_ROWS if r['exclusion'] is not None)),chromosome_J_counts=counts,no_MAF_threshold=True,no_GWAS_outcomes=True,no_native_mask_selection=True))
        rows=[]
        for c in range(1,23):
            pref=root/'qc/controls'/f'chr{c}';rows.append(dict(CHR=c,J_count=counts[str(c)],independent_rows=192,source_window_sha256=sha(str(pref)+'.bed'),B_array_sha256=sha(str(pref)+'.B.npz'),B_receipt_sha256=sha(str(pref)+'.B.receipt.json'),controls=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'],full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
        phase('G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=sha(root/'qc/controls/fixed_G5_selections.json'),chromosomes=rows,maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,finite_reference_503_conditional_only=True,population_LD_or_actual_GWAS_mask_adequacy_established=False,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
        mapping={path:sha(path) for path in controller.expected_worker_outputs(root)}
        common.save(root/'qc/source_worker_result.json',dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=ps,ordinary_output_sha256=mapping,archive_identity=ai,source_index_sha256=sha(root/'extracted/source_index.sqlite'),J_ordered_sha256=QB['logical_sha256']['J_ordered'],chromosome_J_counts=counts,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
        (root/'worker_stdout.log').write_bytes(b'INVENTED_NO_WORKER\n')
        m=dict(status='COMPLETE_REAPED_SOURCE_ONLY',command=command,returncode=0,stop_reason=None,reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256=ps,cleanup_errors=[],run_error=None,process_group_teardown=dict(remaining_group_members=[]),first_teardown_attempt={})
        if fault=='wrong_frozen_command':m['command']=['WRONG_INVENTED_COMMAND']
        if fault=='mutate_invocation_inplace':command.append('--INVENTED_UNAPPROVED_EXTRA_ARG')
        if fault=='extra_directory':(root/'qc/UNREGISTERED_EMPTY_DIRECTORY').mkdir()
        mp=root/'monitor_receipt.json';return mp,common.save(mp,m)
    with patch.object(controller,'exclusive_heavy_lock',no_lock),patch.object(controller,'shared_raw_family_lock',no_lock),patch.object(controller,'identity_gate',lambda *a:None),patch.object(controller,'runtime_profile_gate',lambda *a:None),patch.object(controller,'resource_snapshot',lambda *a,**k:None),patch.object(controller,'Guard',lambda *a:lambda:None),patch.object(controller,'monitor_worker',fake_monitor):
        try:controller.execute(pp,ps,ap,ads);outcome='COMPLETE'
        except RuntimeError as e:outcome='REJECTED: '+str(e)
    committed=not (root/'PENDING.json').exists() and (root/'terminal_seal.json').is_file()
    if fault in ['wrong_frozen_command','extra_directory']:ck(not committed and outcome.startswith('REJECTED') and (root/'PENDING.json').is_file(),'corrected consumer veto')
    else:
        ck(committed and outcome=='COMPLETE','actual metadata terminal complete')
        pr=root/'source_execution_receipt.json';binding=dict(plan_sha256=ps,admission_sha256=ads,scope='SOURCE_AUTHENTICATION_REFERENCE_QC_OPERATOR_ONLY')
        terminal.require_committed(root/'PENDING.json',root/'terminal_seal.json',binding,{str(pr):sha(pr)})
    return dict(outcome=outcome,terminal_committed=committed,initial_controller_argv=initial,final_monitor_argv=json.loads((root/'monitor_receipt.json').read_text())['command'],PENDING_present=(root/'PENDING.json').exists(),actual_worker_or_body_or_mutex_operations=0,normal_producer_qualification='Actual reviewed monitor does not mutate its argv; this is a bounded late-freeze fault witness only.')


def family_meter():
    outer=FIX/'family_meter';outer.mkdir();old=outer/'canonical_signed_reference_source_v1';new=outer/'canonical_signed_reference_source_v2'
    for root in [old,new]:
        root.mkdir();(root/'archive').mkdir();(root/'extracted').mkdir();(root/'qc').mkdir()
    (old/'archive/a').write_bytes(b'123');(old/'extracted/e').write_bytes(b'1234');(old/'qc/q').write_bytes(b'12345');(new/'qc/q').write_bytes(b'12')
    contract=dict(json.loads(PLAN.read_text())['resource_contract']);contract['source_epoch_namespaces']=[str(old),str(new)]
    plan=dict(output_root=str(new),resource_contract=contract,closed_sensitivity_namespace=str(outer/'closed'),compact_representation=compact.SCHEMA,expected_master_reference_rows=1290028)
    got=common.resource_snapshot(plan,time.monotonic())
    ck((got['source_aggregate_bytes'],got['archive_bytes'],got['extraction_bytes'],got['QC_controls_metadata_bytes'])==(14,3,4,7),'old/current family byte accounting')
    return dict(old_and_current_aggregate=14,archive=3,extracted=4,QC=7,only_owned_invented_files_metered=True)


def main():
    ck(shutil.disk_usage('/System/Volumes/Data').free>LIMIT['internal_launch_floor_bytes'],'own internal floor');FIX.mkdir()
    case('actual_bindings_and_inherited_unchanged_science_AST',unchanged)
    case('independent_complete_physical_and_four_logical_reconstruction_binary64',reconstruction)
    case('fixed_capacity_arithmetic_and_preQC_failures',bounds_and_failures)
    case('first_position_q_actual_helper_argv_no_subprocess',curl_command)
    for name,fault in [('healthy',None),('wrong_command','wrong_frozen_command'),('extra_directory','extra_directory'),('inplace_command_rebinding','mutate_invocation_inplace')]:case('actual_controller_'+name,lambda name=name,fault=fault:full_controller(name,fault))
    case('preserved_v1_source_family_in_component_meter',family_meter)
    files={str(p):sha(p) for p in FIX.rglob('*') if p.is_file()}
    result=dict(schema='INDEPENDENT_SIGNED_REFERENCE_SOURCE_V2_DELTA_CONTROLS',actual_plan_sha256=PLAN_SHA,case_count=len(CASES),assertions=ASSERTIONS,cases=CASES,fixture_file_sha256=files,fixture_bytes=sum(Path(p).stat().st_size for p in files),elapsed_seconds=time.monotonic()-START,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,real_source_network_workers_mutex_estimators=0,runtime_census_or_old_suite_runs=0)
    put(R/'independent_signed_reference_source_prelaunch_controls_v2.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['cases','fixture_file_sha256']}))


if __name__=='__main__':main()
