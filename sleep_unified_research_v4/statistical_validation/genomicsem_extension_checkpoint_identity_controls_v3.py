#!/usr/bin/env python3
"""Small synthetic origin/header/family controls; no production bodies read."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
SOURCE=P/'scripts/extension_replay_common_v3.py'
EXPECTED='ca5eaf792eb4a06ef04d323249472452774eada0e54109b645d1b95df6d5bee4'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,sort_keys=True)+'\n');return sha(path)


def family(root):
    root.mkdir();pending=root/'pending.json';seal=root/'seal.json'
    copies=[root/'primary_a.json',root/'primary_b.json']
    paths=[root/('source_%03d.json'%i) for i in range(100)]
    items=[{'path':str(p),'sha256':save(p,{'synthetic_source':i})} for i,p in enumerate(paths)]
    plan={'members':[{'acquisition_receipt':str(p)} for p in paths],
          'acquisition_plan_sha256':'SYNTHETIC_V6_ORIGIN',
          'acquisition_execution_binding':{'executor_sha256':'SYNTHETIC_EXECUTOR'},
          'acquisition_family_terminal_contract':{'pending_path':str(pending),'terminal_seal_path':str(seal),'primary_receipt_paths':[str(p) for p in copies]}}
    primary={'status':'ALL100_EXACT_SOURCE_BODIES_ACQUIRED','plan_sha256':plan['acquisition_plan_sha256'],
             'completed_source_count':100,'reused_exact_source_count':2,'new_exact_source_count':98,
             'source_receipts':items,'pending_path':str(pending),'terminal_seal_path':str(seal)}
    s={'status':'ALL100_FAMILY_TERMINAL_SEAL','plan_sha256':plan['acquisition_plan_sha256'],
       'executor_sha256':'SYNTHETIC_EXECUTOR','source_receipts':items,'pending_path':str(pending),
       'primary_receipt_sha256':{str(p):save(p,primary) for p in copies}}
    save(seal,s)
    return dict(plan=plan,seal=seal,pending=pending,copies=copies,paths=paths,primary=primary,s=s,map={x['path']:x['sha256'] for x in items})


def run():
    assert sha(SOURCE)==EXPECTED
    spec=importlib.util.spec_from_file_location('_checkpoint_identity_v3',SOURCE)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    checks=[]
    def reject(call,code,label):
        try:call()
        except RuntimeError as e:
            assert code in str(e),(label,str(e));checks.append({'case':label,'pass':True,'rejected':code})
        else:raise AssertionError('EXPECTED_REJECTION_MISSING '+label)
    with tempfile.TemporaryDirectory(prefix='checkpoint-identity-') as td:
        tmp=Path(td)
        # Frozen per-origin offsets determine exact command and signed HTTP206
        # content-range admission. These six-byte raw fixtures are synthetic.
        for index in [1,2,3,4]:
            root=tmp/('origin_%d'%index);root.mkdir()
            origin='SYNTHETIC_OLD' if index<=2 else 'SYNTHETIC_V6'
            offset=2 if index in (1,3) else 0
            member={'index':index,'extension_trait_id':'CONTROL_%d'%index,'body_path':str(root/'synthetic_raw'),
                    'expected_bytes':6,'expected_md5':hashlib.md5(b'ABCDEF').hexdigest(),
                    'expected_sha256':hashlib.sha256(b'ABCDEF').hexdigest(),'url':'https://example.invalid/source',
                    'expected_s3_version_id':'SYNTHETIC_VERSION','expected_etag':'SYNTHETIC_ETAG'}
            Path(member['body_path']).write_bytes(b'ABCDEF')
            h={'status':206 if offset else 200,'x-amz-version-id':'SYNTHETIC_VERSION','etag':'SYNTHETIC_ETAG','content-length':str(6-offset)}
            if offset:h['content-range']='bytes 2-5/6'
            header=root/'logs'/(member['extension_trait_id']+'.headers.txt');header.parent.mkdir()
            header.write_text('HTTP/2 '+str(h['status'])+'\n'+''.join(k+': '+v+'\n' for k,v in h.items() if k!='status'))
            identity={'source_folder':str(root),'curl':'/usr/bin/curl','per_body_seconds_limit':100,'resume_offset_by_index':{str(index):offset}}
            plan={'acquisition_operational_identity_by_origin':{origin:identity},'guard':{'internal_floor_bytes':3*1024**3,'SSD_floor_bytes':5*1024**3,'observed_aggregate_worker_RSS_limit_bytes':2*1024**3}}
            receipt=root/'receipt.json';entry={'acquisition_receipt':str(receipt),'acquisition_member':member,'acquisition_origin_plan_sha256':origin}
            logbase=header.with_suffix('').with_suffix('')
            command=['/usr/bin/curl','-q','--fail','--location','--max-redirs','5','--proto','=https','--proto-redir','=https','--tlsv1.2','--max-time','100','--max-filesize','6','--dump-header',str(header),'--output',member['body_path']+'.partial']
            if offset:command+=['--continue-at','2']
            command.append(member['url'])
            state={'internal_free_bytes':3*1024**3,'ssd_free_bytes':5*1024**3}
            r={'status':'EXACT_IMMUTABLE_SOURCE_ACQUIRED','member':member,'plan_sha256':origin,'returncode':0,'stop_reason':None,
               'teardown':{'teardown_verified':True,'remaining_group_members':[]},'post_cleanup_hash_resource_identity_gates_pass':True,
               'actual_size':6,'actual_md5':member['expected_md5'],'actual_sha256':member['expected_sha256'],
               'final_seal_md5':member['expected_md5'],'final_seal_sha256':member['expected_sha256'],'resume_offset':offset,
               'command':command,'observed_headers':h,'resource_before':state,'resource_after':state,
               'headers_sha256':sha(header),'headers_bytes':header.stat().st_size,'peak_observed_owned_rss_bytes':0,'elapsed_seconds':1}
            digest=save(receipt,r);assert m.acquisition_receipt_gate(plan,entry,digest)==r
            checks.append({'case':'exact_origin_command_headers_index%d'%index,'pass':True})
            changed=dict(r,plan_sha256='OTHER_ORIGIN');digest=save(receipt,changed)
            reject(lambda:m.acquisition_receipt_gate(plan,entry,digest),'RAW_SOURCE_NOT_EXACT_SUCCESSFUL_CHECKPOINT','wrong_origin_index%d'%index)
            changed=dict(r,resume_offset=offset+1);digest=save(receipt,changed)
            reject(lambda:m.acquisition_receipt_gate(plan,entry,digest),'ACQUISITION_RESUME_OFFSET_DIFFERS','wrong_offset_index%d'%index)
            if index==3:
                changed=json.loads(json.dumps(r));changed['observed_headers']['status']=200;digest=save(receipt,changed)
                reject(lambda:m.acquisition_receipt_gate(plan,entry,digest),'ACQUISITION_HTTP_IDENTITY_DIFFERS','source3_requires206')
                changed=json.loads(json.dumps(r));changed['observed_headers']['content-range']='bytes 1-5/6';digest=save(receipt,changed)
                reject(lambda:m.acquisition_receipt_gate(plan,entry,digest),'ACQUISITION_EXACT_RANGE_DIFFERS','source3_exact_range_required')
        with patch.object(m,'acquisition_operational_gate',lambda *_:None):
            def case(name,mutate=None,error=None,complete=True,expected_map=True):
                f=family(tmp/('family_%02d'%len(checks)))
                if mutate:mutate(f)
                call=lambda:m.acquisition_family_gate(f['plan'],f['map'] if expected_map else None,require_complete=complete)
                if error:reject(call,error,name)
                else:
                    result=call();assert result is not None
                    checks.append({'case':name,'pass':True})
                return f
            case('exact100_two_identical_primary_copies_and_seal')
            f=family(tmp/'pending_checkpoint');save(f['pending'],{'synthetic_pending':True})
            assert m.acquisition_family_gate(f['plan'],require_complete=False) is None
            checks.append({'case':'individual_checkpoint_credit_allowed_while_family_pending','pass':True})
            reject(lambda:m.acquisition_family_gate(f['plan'],f['map'],require_complete=True),'SOURCE_FAMILY_TERMINAL_PENDING_OR_UNSEALED','pending_vetoes_full_family')
            case('missing_seal',lambda f:f['seal'].unlink(),'SOURCE_FAMILY_TERMINAL_PENDING_OR_UNSEALED')
            case('missing_second_copy',lambda f:f['copies'][1].unlink(),'SOURCE_FAMILY_TWO_EXACT_PRIMARY_COPIES_REQUIRED')
            case('different_second_copy',lambda f:save(f['copies'][1],dict(f['primary'],audit_note='changed')),'SOURCE_FAMILY_TWO_EXACT_PRIMARY_COPIES_REQUIRED')
            def seal_change(f,**updates):f['s'].update(updates);save(f['seal'],f['s'])
            case('wrong_seal_plan',lambda f:seal_change(f,plan_sha256='OTHER'),'SOURCE_FAMILY_TERMINAL_BINDING_DIFFERS')
            case('wrong_seal_executor',lambda f:seal_change(f,executor_sha256='OTHER'),'SOURCE_FAMILY_TERMINAL_BINDING_DIFFERS')
            case('missing_source_entry',lambda f:seal_change(f,source_receipts=f['s']['source_receipts'][:-1]),'SOURCE_FAMILY_FULL100_RECEIPT_SET_DIFFERS')
            case('duplicate_source_entry',lambda f:seal_change(f,source_receipts=f['s']['source_receipts'][:-1]+[f['s']['source_receipts'][0]]),'SOURCE_FAMILY_FULL100_RECEIPT_SET_DIFFERS')
            case('changed_expected_checkpoint_map',lambda f:f['map'].update({next(iter(f['map'])):'OTHER'}),'SOURCE_FAMILY_CHECKPOINT_HASH_MAP_DIFFERS')
            case('changed_source_receipt',lambda f:save(f['paths'][50],{'changed':True}),'SOURCE_FAMILY_FINAL_CHECKPOINT_CHANGED')
            def primary_change(f,**updates):
                f['primary'].update(updates)
                f['s']['primary_receipt_sha256']={str(p):save(p,f['primary']) for p in f['copies']}
                save(f['seal'],f['s'])
            case('wrong_full_count',lambda f:primary_change(f,completed_source_count=99),'SOURCE_FAMILY_PRIMARY_FULL100_BINDING_DIFFERS')
            case('wrong_preserved_count',lambda f:primary_change(f,reused_exact_source_count=1),'SOURCE_FAMILY_PRIMARY_FULL100_BINDING_DIFFERS')
            case('wrong_new_count',lambda f:primary_change(f,new_exact_source_count=97),'SOURCE_FAMILY_PRIMARY_FULL100_BINDING_DIFFERS')
            case('failed_family_primary',lambda f:primary_change(f,status='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'),'SOURCE_ACQUISITION_FAMILY_FAILED_PRESERVED',complete=False)
            case('primary_failure_addendum',lambda f:save(str(f['copies'][0])+'.failure.json',{}),'SOURCE_ACQUISITION_FAMILY_FAILURE_ADDENDUM')
            case('seal_failure_addendum',lambda f:save(str(f['seal'])+'.failure.json',{}),'SOURCE_ACQUISITION_FAMILY_FAILURE_ADDENDUM')
            case('dangling_failure_addendum',lambda f:Path(str(f['seal'])+'.failure.json').symlink_to(f['seal'].parent/'absent'),'SOURCE_ACQUISITION_FAMILY_FAILURE_ADDENDUM')
    assert sha(SOURCE)==EXPECTED
    return {'status':'NARROW_SYNTHETIC_ACQUISITION_ORIGIN_AND_FAMILY_IDENTITY_PASS','source_sha256':EXPECTED,
            'helper_sha256':sha(__file__),'checks':checks,'check_count':len(checks),
            'operational_acquisition_gate_mocked_only_in_family_fixtures':True,
            'real_subprocesses':0,'production_body_reads':0,'real_fits':0,
            'synthetic_body_bytes_per_origin':6,'real_family_completion_claim':False}


if __name__=='__main__':
    out=P/'statistical_validation/genomicsem_extension_checkpoint_identity_controls_receipt_v3.json'
    if out.exists():raise RuntimeError('DISTINCT_CONTROL_RECEIPT_REQUIRED')
    result=run()
    with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','check_count','production_body_reads']}))
