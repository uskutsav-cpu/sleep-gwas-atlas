"""Bounded metadata/guard controls. No real source, worker, network or lock."""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'scripts'));sys.path.insert(0,str(P/'reviews'))
import signed_reference_source_fixture_v1 as fixture
import signed_reference_source_common_v1 as c
import signed_reference_source_controller_v1 as ctrl
import signed_reference_source_worker_v1 as w

RESULTS=[]


def case(name,fn):
    fn();RESULTS.append(dict(case=name,status='PASS'))


def main(out):
    out=Path(out).resolve();out.mkdir();root=out/'invented';root.mkdir();ps='c'*64
    for scenario in ['normal','wrong_version','wrong_UUID','wrong_HEAD','wrong_MD5','redirect']:
        def acquisition(scenario=scenario):
            r=root/('G2_'+scenario);r.mkdir();(r/'qc').mkdir();(r/'archive').mkdir();body=b'INVENTED_NOT_REAL_GENOTYPE_ARCHIVE';md5=hashlib.md5(body).hexdigest()
            cand=dict(file_name='1000G_Phase3_plinkfiles.tgz',file_id='INVENTED_UUID',bytes=len(body),upstream_md5='md5:'+md5,content_url='https://zenodo.org/api/records/8292725/files/1000G_Phase3_plinkfiles.tgz/content')
            p=dict(candidate=cand,curl='INVENTED_NO_CURL_EXECUTION');calls=[]
            old=w.subprocess.run
            def run(command,**kwargs):
                calls.append(command);dest=Path(command[command.index('--output')+1]);url=command[-1]
                assert '--retry' not in command and '--continue-at' not in command and '--location' not in command and '--insecure' not in command
                assert command[command.index('--proto')+1]=='=https';assert command[command.index('--proto-redir')+1]=='=https'
                if url.endswith('/8292725'):
                    value=dict(id=8292725 if scenario!='wrong_version' else 0,doi='10.5281/zenodo.8292725',metadata=dict(access_right='open',license=dict(id='cc-by-4.0')),files=[dict(key=cand['file_name'],id=cand['file_id'] if scenario!='wrong_UUID' else 'WRONG',size=cand['bytes'],checksum=cand['upstream_md5'],links=dict(self=cand['content_url']))]);dest.write_text(json.dumps(value))
                elif '--head' in command:dest.write_text('HTTP/2 200\nContent-Length: '+str(len(body)+(scenario=='wrong_HEAD'))+'\n')
                else:
                    assert command[command.index('--max-filesize')+1]==str(len(body));dest.write_bytes(body if scenario!='wrong_MD5' else b'x'*len(body))
                status='302' if scenario=='redirect' else '200';return SimpleNamespace(returncode=0,stdout=(status+'\n'+url+'\n').encode(),stderr=b'')
            class Guard:
                def __call__(self,**kw):pass
            w.subprocess.run=run
            try:
                if scenario=='normal':archive,identity=w.acquire(p,r,Guard());assert identity['md5']==md5;assert len(calls)==3
                else:fixture.reject(lambda:w.acquire(p,r,Guard()))
            finally:w.subprocess.run=old
        case('mock_exact_source_G2_'+scenario,acquisition)
    for field in ['missing_output','symlink_same_bytes','rebound_result','G3_samples','G4_rank_count','G5_missing_chromosome','G5_B_padding_claim','G5_C_credit']:
        def consumption(field=field):
            r=root/field;p,result,m=fixture.fake_source_metadata(r,ps);mp=r/'monitor.json';mh=c.save(mp,m);rp=r/'qc/source_worker_result.json';rh=c.regular_sha(rp)
            if field=='missing_output':Path(next(iter(result['ordinary_output_sha256']))).unlink()
            elif field=='symlink_same_bytes':
                file=r/'qc/static_selection.jsonl';payload=file.read_bytes();outside=root/(field+'.outside');outside.write_bytes(payload);file.unlink();file.symlink_to(outside)
            elif field=='rebound_result':
                result['additional']='INVENTED_REBIND';rp.write_text(json.dumps(result,indent=2)+'\n')
            else:
                name='G3_format.json' if field=='G3_samples' else 'G4_static_seal.json' if field=='G4_rank_count' else 'G5_operator_controls.json';fp=r/'qc'/name;value=json.loads(fp.read_text())
                if field=='G3_samples':value['samples']=502
                elif field=='G4_rank_count':value['master_reference_rows']-=1
                elif field=='G5_missing_chromosome':value['chromosomes'].pop()
                elif field=='G5_B_padding_claim':value['maintained_B_padding_independent_rejection']=True
                elif field=='G5_C_credit':value['C_calibration_credit']=1
                fp.write_text(json.dumps(value,indent=2)+'\n');result['ordinary_output_sha256'][str(fp)]=c.regular_sha(fp);rp.write_text(json.dumps(result,indent=2)+'\n');rh=c.regular_sha(rp)
            fixture.reject(lambda:ctrl.consume_worker_result(p,mp,mh,ps,rh))
        case('strict_consumption_'+field,consumption)
    def meter_controls():
        p=dict(output_root=str(root/'meter'),closed_sensitivity_namespace=str(root/'closed'),resource_contract=dict(internal_floor_bytes=3<<30,SSD_launch_free_required_bytes=8<<30,SSD_floor_bytes=5<<30,source_aggregate_cap_bytes=3<<30,archive_and_all_failed_partial_cap_bytes=512<<20,extracted_all_partial_and_failed_cap_bytes=2<<30,QC_controls_metadata_and_all_failed_output_cap_bytes=512<<20,global_ceiling_bytes=300<<30,closed_sensitivity_cap_bytes=256<<20,whole_stage_monotonic_deadline_seconds=7200))
        oldtree=c.full_tree_bytes;olddisk=c.shutil.disk_usage
        try:
            count={p['output_root']:0,p['output_root']+'/archive':0,p['output_root']+'/extracted':0,str(Path(p['output_root']).parent):0,p['closed_sensitivity_namespace']:0}
            c.full_tree_bytes=lambda path:count[str(path)];c.shutil.disk_usage=lambda path:SimpleNamespace(free=100<<30)
            c.resource_snapshot(p,time.monotonic(),True)
            for key,value in [(p['output_root'],(3<<30)+1),(p['output_root']+'/archive',(512<<20)+1),(p['output_root']+'/extracted',(2<<30)+1),(str(Path(p['output_root']).parent),(300<<30)+1),(p['closed_sensitivity_namespace'],(256<<20)+1)]:
                count[key]=value;fixture.reject(lambda:c.resource_snapshot(p,time.monotonic()));count[key]=0
            c.shutil.disk_usage=lambda path:SimpleNamespace(free=(3<<30)-1);fixture.reject(lambda:c.resource_snapshot(p,time.monotonic()))
            c.shutil.disk_usage=lambda path:SimpleNamespace(free=4<<30 if str(path)!='/System/Volumes/Data' else 100<<30);fixture.reject(lambda:c.resource_snapshot(p,time.monotonic()))
            c.shutil.disk_usage=lambda path:SimpleNamespace(free=7<<30 if str(path)!='/System/Volumes/Data' else 100<<30);fixture.reject(lambda:c.resource_snapshot(p,time.monotonic(),True))
            c.shutil.disk_usage=lambda path:SimpleNamespace(free=100<<30);fixture.reject(lambda:c.resource_snapshot(p,time.monotonic()-7201))
        finally:c.full_tree_bytes=oldtree;c.shutil.disk_usage=olddisk
    case('all_source_subcaps_global_closed_floors_launch_and_deadline',meter_controls)
    def raw_lifetime():
        p=root/'FAKE_RAW_LOCK_FILE_NOT_REAL_MUTEX';p.write_bytes(b'INVENTED');events=[];oldflock=c.fcntl.flock;oldsleep=c.time.sleep
        c.fcntl.flock=lambda fd,operation:events.append(operation);c.time.sleep=lambda _:None
        attempts=[]
        def clear():
            attempts.append(1)
            if len(attempts)==1:raise RuntimeError('INVENTED_UNVERIFIED_CLEANUP')
        try:
            with c.shared_raw_family_lock(p,clear):assert events==[c.fcntl.LOCK_EX|c.fcntl.LOCK_NB]
            assert len(attempts)==2;assert events==[c.fcntl.LOCK_EX|c.fcntl.LOCK_NB,c.fcntl.LOCK_UN]
        finally:c.fcntl.flock=oldflock;c.time.sleep=oldsleep
    case('mock_raw_mutex_release_only_after_verified_retry',raw_lifetime)
    def cleanup_failure():
        oldhelper=ctrl.monitor_helper;oldgroup=ctrl.group_exists;oldkill=ctrl.os.killpg;events=[]
        proc=SimpleNamespace(pid=99999999,wait=lambda **kw:0);entry=dict(proc=proc,cleanup_errors=[]);owned=[entry]
        ctrl.monitor_helper=lambda plan:SimpleNamespace(terminate_owned=lambda proc:(_ for _ in ()).throw(RuntimeError('INVENTED_HELPER_FAILURE')))
        states=iter([True,False]);ctrl.group_exists=lambda pid:next(states);ctrl.os.killpg=lambda pid,sign:events.append((pid,sign))
        try:ctrl.cleanup({},owned);assert not owned and entry['cleanup_errors'];assert events==[(99999999,ctrl.signal.SIGKILL)]
        finally:ctrl.monitor_helper=oldhelper;ctrl.group_exists=oldgroup;ctrl.os.killpg=oldkill
    case('mock_owned_cleanup_independent_fallback_retains_first_error',cleanup_failure)
    plan=P/'manifests/signed_reference_source_execution_plan_v1.json';v=c.load_frozen(plan,'6d73f45ee3f73f7db7c0cc727d9b2d3e1af30d3ee27391cf1154f1f8ab8482b4')
    for path,digest in v['dependency_sha256'].items():
        if path not in v['preparation_deferred_dependency_sha256']:assert c.regular_sha(path)==digest
    receipt=dict(status='AUTHOR_METADATA_SOURCE_CONTROLS_PASS_NO_EXECUTION_ADMISSION',case_count=len(RESULTS),cases=RESULTS,
        frozen_plan_sha256=c.regular_sha(plan),fresh_nonbody_dependency_count=len(v['dependency_sha256'])-len(v['preparation_deferred_dependency_sha256']),
        inherited_deferred_dependency_count=26,actual_source_body_reads=0,actual_network_workers_mutexes=0,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,
        mock_completion_is_metadata_only_not_G2_G5_achievement=True)
    c.save(out/'metadata_control_receipt.json',receipt);print(json.dumps(dict(status=receipt['status'],cases=len(RESULTS))))


if __name__=='__main__':main(sys.argv[1])
