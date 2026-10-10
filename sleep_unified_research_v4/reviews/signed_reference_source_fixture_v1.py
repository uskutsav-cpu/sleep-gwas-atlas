"""Invented source-only fixtures. Does not open a real BED/BIM/FAM/archive."""
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tarfile
import tempfile
import time
from contextlib import contextmanager

S=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(S))
import numpy as np
import signed_reference_decode_v1 as a
import signed_reference_source_common_v1 as common
import signed_reference_source_controller_v1 as controller

PYB='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python'
READER='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/ldsc/ldscore/ldscore.py'
CASES=[]


def check(name,fn):
    try:fn();CASES.append(dict(case=name,status='PASS'))
    except BaseException as error:CASES.append(dict(case=name,status='FAIL',error=repr(error)));raise


def reject(fn,contains=None):
    try:fn()
    except BaseException as e:
        if contains is not None:assert contains in str(e),repr(e)
        return
    raise AssertionError('Expected rejection')


def packed(dose):
    values=np.array([0 if x==2 else 1 if x==-1 else 2 if x==1 else 3 for x in dose]+[0],dtype=np.uint8)
    return bytes(np.sum(values.reshape(126,4)*np.array([1,4,16,64],dtype=np.uint8),axis=1).astype(np.uint8))


def archive_fixture(path,change=None,prefix=''):
    buf=io.BytesIO()
    with tarfile.open(fileobj=buf,mode='w',format=tarfile.USTAR_FORMAT) as tar:
        for c in range(1,23):
            for ext in ['bed','bim','fam']:
                name=prefix+f'1000G.EUR.QC.{c}.{ext}';item=tarfile.TarInfo(name);body=b'invented\n';item.size=len(body)
                if change and (c,ext)==(1,'bed'):item,body=change(item,body)
                tar.addfile(item,io.BytesIO(body) if body else None)
    path.write_bytes(gzip.compress(buf.getvalue(),mtime=0));return buf.getvalue()


def fake_source_metadata(root,ps):
    """Metadata stubs expressly earn no G2--G5 science/body credit."""
    root=Path(root);root.mkdir(exist_ok=True);(root/'worker_stdout.log').write_bytes(b'')
    mapping={}
    for p in sorted(controller.expected_worker_outputs(root)):
        p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'INVENTED_METADATA_STUB\n');mapping[str(p)]=common.regular_sha(p)
    candidate=dict(bytes=21,upstream_md5='md5:INVENTED',file_name='INVENTED')
    archive=root/'archive/1000G_Phase3_plinkfiles.tgz';identity=dict(bytes=21,md5='INVENTED',sha256=mapping[str(archive)])
    counts={str(c):211 for c in range(1,23)};plan=dict(output_root=str(root),candidate=candidate,expected_master_reference_rows=1290028)
    def put(name,value):
        p=root/'qc'/name;p.write_text(json.dumps(value,indent=2)+'\n');mapping[str(p)]=common.regular_sha(p)
    put('G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',candidate=candidate,archive_identity=identity))
    put('G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,chromosome_variant_counts=counts,source_index_sha256=mapping[str(root/'extracted/source_index.sqlite')]))
    put('G4_static_seal.json',dict(status='FROZEN_BEFORE_GENOTYPE_QC',master_reference_rows=1290028,static_ordered_sha256=mapping[str(root/'qc/static_selection.jsonl')],static_exclusions_sha256=mapping[str(root/'qc/static_exclusions.jsonl')]))
    put('G4_J_seal.json',dict(status='REFERENCE_ONLY_J_FROZEN',J_ordered_sha256=mapping[str(root/'qc/J_ordered.jsonl')],chromosome_J_counts=counts,J_count=sum(counts.values()),no_MAF_threshold=True,no_GWAS_outcomes=True,no_native_mask_selection=True))
    members=[]
    for c in range(1,23):
        pref=root/'qc/controls'/f'chr{c}';members.append(dict(CHR=c,J_count=211,independent_rows=211,source_window_sha256=mapping[str(pref)+'.bed'],B_array_sha256=mapping[str(pref)+'.B.npz'],B_receipt_sha256=mapping[str(pref)+'.B.receipt.json'],controls=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'],full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
    put('G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=mapping[str(root/'qc/controls/fixed_G5_selections.json')],chromosomes=members,maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,finite_reference_503_conditional_only=True,population_LD_or_actual_GWAS_mask_adequacy_established=False,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
    result=dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=ps,ordinary_output_sha256=mapping,archive_identity=identity,source_index_sha256=mapping[str(root/'extracted/source_index.sqlite')],J_ordered_sha256=mapping[str(root/'qc/J_ordered.jsonl')],chromosome_J_counts=counts,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False)
    p=root/'qc/source_worker_result.json';common.save(p,result)
    monitor=dict(status='COMPLETE_REAPED_SOURCE_ONLY',returncode=0,stop_reason=None,reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256=ps,cleanup_errors=[],run_error=None,process_group_teardown=dict(remaining_group_members=[]),first_teardown_attempt={})
    return plan,result,monitor


def main(out):
    if sys.version_info[:3]!=(3,11,11) or np.__version__!='1.26.4':raise RuntimeError('DECLARED_FIXTURE_RUNTIME_REQUIRED')
    out=Path(out).resolve();out.mkdir();start=time.monotonic();tmp=out/'invented_private';tmp.mkdir()
    dose=np.array(([2,-1,1,0]*126)[:503],dtype=np.int8)
    good=np.array(([2,0,1,2,1,0]*84)[:503],dtype=np.int8)
    def exact_states():assert np.array_equal(a.decode_a1(packed(dose),np)[:,0],dose)
    check('A_exact503_all_four_states',exact_states)
    check('A_nonzero_padding_reject',lambda:reject(lambda:a.decode_a1(packed(good)[:-1]+bytes([packed(good)[-1]|0xC0]),np),'PADDING'))
    check('A_partial_row_reject',lambda:reject(lambda:a.decode_a1(packed(good)[:-1],np),'COMPLETE'))
    check('A_empty_chunk_reject',lambda:reject(lambda:a.decode_a1(b'',np),'CHUNK'))
    check('A_gt2048_chunk_reject',lambda:reject(lambda:a.decode_a1(packed(good)*2049,np),'CHUNK'))
    def boundaries():
        x=good.copy();x[:25]=-1;assert a.qc_column(x,np)['exclusion'] is None
        x[25]=-1;assert a.qc_column(x,np)['exclusion']=='MISSING_GT25_OF503'
        for x,reason in [(np.full(503,-1),'ALL_MISSING'),(np.r_[0,np.full(502,-1)],'FEWER_THAN_TWO_OBSERVED'),(np.ones(503),'EXACT_ZERO_VARIANCE')]:assert a.qc_column(x,np)['exclusion']==reason
        x=np.zeros(503,dtype=np.int8);x[-1]=1;assert a.qc_column(x,np)['integer_variance_numerator']==502;assert a.qc_column(x,np)['exclusion'] is None
    check('A_missingness_exact_integer_variance_boundaries',boundaries)
    def normalization():
        x=good.copy();x[:25]=-1;y=a.standardized(x,1,np);assert abs(y.mean())<1e-10;assert abs(y@y/503-1)<1e-10;assert np.array_equal(-y,a.standardized(x,-1,np))
        reject(lambda:a.standardized(np.ones(503),1,np),'ZERO_VARIANCE')
        assert a.allele_sign(('A','C'),('A','C'))==(1,None)
        assert a.allele_sign(('A','C'),('C','A'))==(-1,None)
        assert a.allele_sign(('A','C'),('T','G'))[0] is None
        assert a.allele_sign(('A','T'),('A','T'))[1]=='PALINDROMIC'
    check('A_mean_imputation_signed_orientation_no_floor',normalization)
    def selections():
        reject(lambda:a.frozen_control_indices(191),'INCOMPLETE')
        for n in [192,193,288,10001]:
            v=a.frozen_control_indices(n);assert len(set(v['first64']+v['middle64']+v['last64']))==192;assert len(set(v['early32_scatter']+v['middle32_scatter']+v['late32_scatter']))==96
    check('fixed_windows_scatter_incomplete_and_distinct',selections)
    def operators():
        G=np.column_stack([a.standardized(np.roll(good,i),(-1)**i,np) for i in range(211)])
        rows=list(range(211));v=np.array([(i%13-6)/6 for i in rows]);w=np.array([(i%17-8)/8 for i in rows]);R=G.T@G/503
        read=lambda ix:G[:,ix]
        for c in [1,7,64,2048]:
            for rev in [False,True]:
                got,gv=a.operator(rows,read,v,np,c,rev);assert a.close_array(got,R@v,np);assert a.close_array(np.array([v@got]),np.array([gv@gv/503]),np)
        rw,_=a.operator(rows,read,w,np);assert a.close_array(np.array([v@rw]),np.array([w@(R@v)]),np)
        d=np.where(np.arange(211)%3==0,-1.,1.);assert a.close_array((G*d).T@(G*d)/503,d[:,None]*R*d[None,:],np)
        broken=R.copy();broken[:64,64:]=0;broken[64:,:64]=0;assert not a.close_array(broken@v,R@v,np)
        assert np.linalg.matrix_rank(G)<=502
    check('signed_Gram_independent_arithmetic_chunks_reverse_distant_operator',operators)
    def normal_tar():
        for pref in ['', '1000G_EUR_Phase3_plink/']:
            p=tmp/f'normal{len(pref)}.tgz';archive_fixture(p,prefix=pref);r=a.extract_strict66(p,tmp/f'extract{len(pref)}',100000);assert len(r['files'])==66
    check('manual_strict66_both_registered_prefixes',normal_tar)
    def kind(k,name=None):
        def change(item,body):
            item.type=k
            if name:item.name=name
            if k!=tarfile.REGTYPE:item.size=0;item.linkname='target';body=None
            return item,body
        return change
    for label,change in [('traversal',kind(tarfile.REGTYPE,'../escape')),('absolute',kind(tarfile.REGTYPE,'/escape')),('symlink',kind(tarfile.SYMTYPE)),('hardlink',kind(tarfile.LNKTYPE)),('fifo',kind(tarfile.FIFOTYPE)),('unexpected',kind(tarfile.REGTYPE,'unexpected')),('duplicate',kind(tarfile.REGTYPE,'1000G.EUR.QC.1.bim')),('mixedprefix',kind(tarfile.REGTYPE,'1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bed'))]:
        def case(change=change,label=label):
            p=tmp/(label+'.tgz');archive_fixture(p,change);reject(lambda:a.extract_strict66(p,tmp/(label+'extract'),100000))
        check('archive_'+label+'_reject',case)
    def gzip_faults():
        base=tmp/'gzbase.tgz';archive_fixture(base)
        for label,data in [('truncated',base.read_bytes()[:-4]),('CRC',base.read_bytes()[:-8]+b'\0'*8),('trailing',base.read_bytes()+b'extra'),('two_members',base.read_bytes()+gzip.compress(b''))]:
            p=tmp/(label+'.tgz');p.write_bytes(data);reject(lambda:a.extract_strict66(p,tmp/(label+'extract'),100000))
        reject(lambda:a.extract_strict66(base,tmp/'capextract',1),'CAP')
    check('single_gzip_CRC_EOF_trailing_concat_and_extraction_cap',gzip_faults)
    def maintained_B():
        cols=[good,np.roll(good,13),np.roll(good,29),np.where(np.arange(503)<25,-1,np.roll(good,43))]
        bed=tmp/'B_good.bed';bed.write_bytes(a.BED_HEADER+b''.join(packed(x) for x in cols));rows=[dict(SNP=f'INVENTED{i}',CHR=1,BP=10+i) for i in range(4)]
        spec=dict(samples=503,rows=rows,bed_path=str(bed),bed_sha256=hashlib.sha256(bed.read_bytes()).hexdigest(),maintained_reader_path=READER,output_path=str(tmp/'Bgood.npz'),receipt_path=str(tmp/'Bgood.json'))
        p=tmp/'B_metadata.json';h=common.save(p,spec)
        env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
        cmd=[PYB,'-B',str(S/'signed_reference_ldsc_decoderB_v1.py'),'--metadata',str(p),'--metadata-sha256',h]
        result=subprocess.run(cmd,capture_output=True,text=True,env=env);(tmp/'Bgood.log').write_text(result.stdout+result.stderr);assert result.returncode==0,result.stderr
        receipt=common.load_frozen(spec['receipt_path'],common.regular_sha(spec['receipt_path']))
        with np.load(spec['output_path'],allow_pickle=False) as data:raw=data['A2_dose'];Y=data['standardized_A2']
        assert np.array_equal(np.column_stack(cols),np.where(raw==9,-1,2-raw))
        G=np.column_stack([a.standardized(x,1,np) for x in cols]);assert a.close_array(G,-Y,np)
        assert a.close_array(G*np.array([1,-1,1,-1]),-Y*np.array([1,-1,1,-1]),np)
        assert receipt['maintained_reader_sha256']=='df98f324ca775ece05a6e8f93e0584c88bde9535c360fda6ff96bb7f475e726b'
        for label,payload,want in [('padding',a.BED_HEADER+packed(good)[:-1]+bytes([packed(good)[-1]|0xC0]),0),('zero',a.BED_HEADER+packed(np.ones(503)),1),('allmissing',a.BED_HEADER+packed(np.full(503,-1)),1),('samplemajor',b'\x6c\x1b\x00'+packed(good),1),('truncated',a.BED_HEADER+packed(good)[:-1],1),('trailing',a.BED_HEADER+packed(good)+b'x',1)]:
            b=tmp/(label+'B.bed');b.write_bytes(payload);q=tmp/(label+'B.metadata.json');sp={**spec,'rows':rows[:1],'bed_path':str(b),'bed_sha256':hashlib.sha256(payload).hexdigest(),'output_path':str(tmp/(label+'B.npz')),'receipt_path':str(tmp/(label+'B.json'))};qh=common.save(q,sp)
            rr=subprocess.run([PYB,'-B',str(S/'signed_reference_ldsc_decoderB_v1.py'),'--metadata',str(q),'--metadata-sha256',qh],capture_output=True,text=True,env=env);(tmp/(label+'B.log')).write_text(rr.stdout+rr.stderr)
            assert (rr.returncode==0)==(want==0),(label,rr.stderr)
        # Padding acceptance is recorded as a library qualification, never PASS
        # for G3. Original zero/filter rejections give no alternate J rule.
    check('maintained_B_actual_original_loader_orientation_missingness_filter_padding_limits',maintained_B)
    def consumers():
        root=tmp/'consumer';ps='a'*64;plan,result,goodm=fake_source_metadata(root,ps)
        for name,field,value in [('healthy',None,None),('cleanup','cleanup_errors',['bad']),('run','run_error',{}),('teardown','process_group_teardown',dict(remaining_group_members=[],cleanup_error='bad')),('first','first_teardown_attempt',dict(exception=None))]:
            m=root/(name+'.monitor.json');h=common.save(m,{**goodm,**({field:value} if field else {})})
            if field:reject(lambda:controller.consume_worker_result(plan,m,h,ps))
            else:controller.consume_worker_result(plan,m,h,ps)
        ordinary=root/'qc/static_selection.jsonl';ordinary.write_bytes(b'mutated');m=root/'healthy.monitor.json';reject(lambda:controller.consume_worker_result(plan,m,common.regular_sha(m),ps),'CHANGED')
    check('source_monitor_all_four_failure_fields_and_output_mutation',consumers)
    def namespace_companions():
        root=tmp/'registered';root.mkdir();p=root/'data';p.write_bytes(b'INVENTED')
        ad=root/'._data';ad.write_bytes(struct.pack('>II16sH',0x00051607,0x00020000,b'\0'*16,1)+struct.pack('>III',2,38,1)+b'x')
        frozen={str(p):common.regular_sha(p)};common.freeze_registered(root,[root,p],frozen);assert str(ad) in frozen
        (root/'unregistered').write_bytes(b'x');reject(lambda:common.freeze_registered(root,[root,p],frozen),'INVENTORY')
    check('exact_registered_AppleDouble_and_unexpected_extra_reject',namespace_companions)
    def source_format_controls():
        import signed_reference_source_worker_v1 as worker
        root=tmp/'format';root.mkdir();(root/'qc').mkdir();(root/'extracted').mkdir();members=root/'extracted/members';members.mkdir()
        sample=[f'INVENTED_SAMPLE{i:03}' for i in range(503)];fam=''.join(f'INVENTED_FID {s} 0 0 0 -9\n' for s in sample).encode()
        panel=root/'panel.tsv';panel.write_text('sample\tpop\tsuper_pop\tgender\n'+''.join(f'{s}\tINVENTED\tEUR\tINVENTED\n' for s in sample));files={}
        for c in range(1,23):
            for ext,payload in [('fam',fam),('bim',f'{c} INVENTED_CHR{c} 0 100 A C\n'.encode()),('bed',a.BED_HEADER+packed(good))]:
                p=members/f'1000G.EUR.QC.{c}.{ext}';p.write_bytes(payload);files[p.name]=dict(path=str(p),sha256=common.regular_sha(p))
        plan=dict(authoritative_panel_path=str(panel),reference_body_sha256={str(panel):common.regular_sha(panel)})
        class Guard:
            def __call__(self,**kw):pass
        source,sh,paths=worker.source_format(plan,root,dict(files=files),Guard());assert len(paths)==22;assert common.regular_sha(source)==sh
        coordinate=root/'coordinate.tsv.gz';coord='SNP\tCHR\tBP\tA1\tA2\n';references=[]
        for c in range(1,23):
            SNP=f'INVENTED_CHR{c}' if c!=5 else 'INVENTED_MISSING_SOURCE';bp=25000000 if c==6 else 100
            r=root/f'{c}.l2.ldscore.gz';r.write_bytes(gzip.compress(f'CHR SNP BP L2\n{c} {SNP} {bp} 2\n'.encode(),mtime=0));references.append(str(r));plan['reference_body_sha256'][str(r)]=common.regular_sha(r)
            if c!=4:coord+=f'{SNP}\t{c}\t{101 if c==2 else bp}\t'+('C\tA' if c==1 else 'G\tA' if c==3 else 'A\tC')+'\n'
        coordinate.write_bytes(gzip.compress(coord.encode(),mtime=0));plan['coordinate_map_path']=str(coordinate);plan['reference_body_sha256'][str(coordinate)]=common.regular_sha(coordinate);plan['reference_files_numeric_order']=references;plan['expected_master_reference_rows']=22
        static,ss=worker.prepare_static(plan,root,source,sh,Guard());selected=list(worker.consume_jsonl(static,ss,Guard()))
        assert len(selected)==17;assert selected[0]['master_rank0']==0 and selected[0]['sign']==-1
        assert selected[1]['master_rank0']==6 and selected[1]['master_rank1']==7
        ex=[json.loads(x)['exclusion'] for x in (root/'qc/static_exclusions.jsonl').read_text().splitlines()]
        assert ex==['REFERENCE_COORDINATE_MAP_CONFLICT','ALLELE_CONFLICT_NO_COMPLEMENT_RESCUE','MISSING_COORDINATE_MAP_ID','MISSING_SOURCE_ID','MHC_INCLUSIVE']
        reject(lambda:worker.genotype_qc(root,static,ss,paths,np,Guard()),'FEWER_THAN192')
        j=[json.loads(x) for x in (root/'qc/J_ordered.jsonl').read_text().splitlines()]
        assert len(j)==17;assert j[1]['master_rank0']==6 and j[1]['J_index0']==1
        # Separate fresh roots are needed; failed files/receipts are retained.
        for label,ext,payload in [('wrong_fam_order','fam',b'\n'.join(fam.splitlines()[::-1])+b'\n'),('badpadding','bed',a.BED_HEADER+packed(good)[:-1]+bytes([packed(good)[-1]|0xC0])),('duplicate_source_id','bim',b'2 INVENTED_CHR1 0 100 A C\n')]:
            r=tmp/('format_'+label);r.mkdir();(r/'qc').mkdir();(r/'extracted').mkdir();changed={k:dict(v) for k,v in files.items()};p=r/f'1000G.EUR.QC.2.{ext}';p.write_bytes(payload);changed[p.name]=dict(path=str(p),sha256=common.regular_sha(p));reject(lambda:worker.source_format(plan,r,dict(files=changed),Guard()))
    check('all22_503_sample_order_BIM_global_unique_and_BED_padding_format',source_format_controls)
    def runtime_link_contracts():
        envA=tmp/'runtimeA';envB=tmp/'runtimeB'
        def profile(env,strict):
            env.mkdir();(env/'real').mkdir();f=env/'real/data';f.write_bytes(b'INVENTED_RUNTIME');(env/'link').symlink_to('real',target_is_directory=True)
            regulars={str(f):dict(bytes=f.stat().st_size,sha256=common.regular_sha(f))}
            dirs=sorted([str(env),str(env/'real')]);tree=dict(regular_files=regulars,physical_directories=[str(env/'real')],file_symlinks={},directory_symlinks={})
            item=dict(literal_link='real',resolved_target=str(env/'real'))
            p=dict(environment_root=str(env),regular_files=regulars,symlinks={})
            if strict:p.update(physical_directories=dirs,symlink_directories={str(env/'link'):{**item,'tree_sha256':hashlib.sha256(json.dumps(tree,sort_keys=True,separators=(',',':')).encode()).hexdigest()}})
            else:p['symlinks'][str(env/'link')]={**item,'target_is_file':False}
            out=tmp/(env.name+'.profile.json');h=common.save(out,p);return out,h,p
        ap,ah,A=profile(envA,False);bp,bh,B=profile(envB,True)
        p=dict(runtime_A_profile_path=str(ap),runtime_A_profile_sha256=ah,runtime_A_environment_root=str(envA),runtime_B_required_profile_path=str(bp),runtime_B_profile_sha256=bh,runtime_B_environment_root=str(envB))
        admission=dict(runtime_B_profile_sha256=bh,independent_complete_B_runtime_profile_review_pass=True)
        controller.runtime_profile_gate(p,admission)
        reject(lambda:controller.runtime_profile_gate(p,{**admission,'independent_complete_B_runtime_profile_review_pass':False}),'PENDING')
        (envB/'unregistered').mkdir();reject(lambda:controller.runtime_profile_gate(p,admission),'DIRECTORY_CENSUS');(envB/'unregistered').rmdir()
        (envB/'real/data').write_bytes(b'CORRUPTED');reject(lambda:controller.runtime_profile_gate(p,admission),'CHANGED')
    check('runtime_B_contained_directory_tree_and_A_legacy_link_qualification',runtime_link_contracts)
    def terminal_stubs():
        # Patch both mutex entrypoints and all worker/runtime/physical/resource
        # operations. This tests only the real metadata controller/Terminal2.
        original={name:getattr(controller,name) for name in ['identity_gate','runtime_profile_gate','resource_snapshot','exclusive_heavy_lock','shared_raw_family_lock','monitor_worker','cleanup']}
        @contextmanager
        def lock(*args,**kw):yield 17
        def cleanup(*args):assert not args[1]
        try:
            controller.identity_gate=lambda *a:None;controller.runtime_profile_gate=lambda *a:None
            controller.exclusive_heavy_lock=lock;controller.shared_raw_family_lock=lock;controller.cleanup=cleanup
            for label in ['healthy','postseal_resource','changed_after_parse']:
                root=tmp/('terminal_'+label);p=tmp/(label+'.plan.json');ad=tmp/(label+'.admit.json');heavy=tmp/'fake_heavy';heavy.write_bytes(b'INVENTED_NOT_A_MUTEX')
                plan=dict(output_root=str(root),resource_contract=dict(shared_heavy_mutex=str(heavy)),raw_transfer_family_mutex=str(heavy),source_python='INVENTED_NO_PROCESS',source_worker='INVENTED_NO_WORKER',candidate=dict(bytes=21,upstream_md5='md5:INVENTED',file_name='INVENTED'),expected_master_reference_rows=1290028)
                ph=common.save(p,plan);ah=common.save(ad,{})
                def monitor(plan,command,fds,started,ownership,gates):
                    _,_,m=fake_source_metadata(root,ph);mfile=root/'monitor_receipt.json';mh=common.save(mfile,m);return mfile,mh
                controller.monitor_worker=monitor
                def resource(plan,started,launch=False):
                    if label=='postseal_resource' and (root/'terminal_seal.json').exists():raise RuntimeError('INVENTED_FINAL_RESOURCE_FLOOR')
                    return {}
                controller.resource_snapshot=resource
                oldload=controller.load_frozen
                if label=='changed_after_parse':
                    def changed(path,h):
                        value=oldload(path,h)
                        if Path(path)==root/'qc/G3_format.json' and (root/'source_execution_receipt.json').exists():
                            (root/'qc/static_selection.jsonl').write_bytes(b'INVENTED_POSTPARSE_MUTATION')
                        return value
                    controller.load_frozen=changed
                try:
                    if label=='healthy':controller.execute(p,ph,ad,ah);assert not (root/'PENDING.json').exists();assert (root/'terminal_seal.json').exists()
                    else:reject(lambda:controller.execute(p,ph,ad,ah));assert (root/'PENDING.json').exists();assert (root/'source_execution_receipt.json.failure.json').exists()
                finally:controller.load_frozen=oldload
        finally:
            for name,value in original.items():setattr(controller,name,value)
    check('whole_controller_metadata_stub_Terminal2_postpersist_resource_and_parse_mutation',terminal_stubs)
    result=dict(status='AUTHOR_INVENTED_FIXTURES_PASS_NO_PRODUCTION_ADMISSION',cases=CASES,case_count=len(CASES),
        elapsed_seconds=time.monotonic()-start,python=sys.version,numpy=np.__version__,
        real_genotype_or_archive_body_reads=0,real_worker_or_mutex_or_network_operations=0,
        decoderB_original_padding_rejection_not_provided=True,full_B_runtime_audit_pending=True,
        S_simulations_run=0,C_calibration_credit=0,estimator_calls=0)
    (out/'fixture_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],cases=len(CASES))))


if __name__=='__main__':main(sys.argv[1])
