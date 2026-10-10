"""Invented compact equivalence/curl/consumer delta; no actual source operation."""
import ast
import contextlib
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import sqlite3
import shutil
import struct
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'scripts'))
import signed_reference_source_worker_v1 as old
import signed_reference_source_worker_v2 as new
import signed_reference_source_common_v1 as c1
import signed_reference_source_common_v2 as c2
import signed_reference_source_controller_v2 as controller
import signed_reference_compact_v2 as compact
import numpy as np

OUT=Path(sys.argv[1])
cases=[]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    b=(json.dumps(v,sort_keys=True,indent=2)+'\n').encode();h=hashlib.sha256(b).hexdigest()
    with Path(p).open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
    assert sha(p)==h;return h
def expect(name,fn):
    try:fn()
    except (RuntimeError,ValueError,sqlite3.IntegrityError):cases.append(dict(case=name,status='PASS'));return
    raise AssertionError('EXPECTED_REJECT '+name)
def pass_case(name):cases.append(dict(case=name,status='PASS'))
def original_rows(path):return [json.loads(s) for s in Path(path).read_text().splitlines()]
def encode(dose):
    values={2:0,-1:1,1:2,0:3};out=bytearray(126)
    for i,n in enumerate(dose):out[i//4]|=values[int(n)]<<((i%4)*2)
    return bytes(out)
def dirs(root):
    root.mkdir()
    for d in ['archive','extracted','extracted/members','qc','qc/controls','tmp']:(root/d).mkdir()


def main():
    assert sys.version_info[:3]==(3,11,11) and np.__version__=='1.26.4'
    OUT.mkdir();assets=OUT/'INVENTED_ONLY';assets.mkdir()
    refs=assets/'refs';refs.mkdir();coordinate=assets/'coordinates.tsv.gz'
    source=assets/'source.sqlite';db=sqlite3.connect(source)
    db.execute('CREATE TABLE source(SNP TEXT PRIMARY KEY,CHR INT,BP INT,A1 TEXT,A2 TEXT,row0 INT,offset INT)')
    reference_paths=[];coordinates=[];beds={};total=0;refhash={}
    base=np.array([i%3 for i in range(503)],dtype=np.int8)
    genotype_special=[np.full(503,-1,dtype=np.int8),np.r_[np.array([1]),np.full(502,-1)],
        np.zeros(503,dtype=np.int8),np.r_[np.full(26,-1),base[26:]],
        np.r_[np.full(25,-1),base[25:]],np.r_[np.array([1]),np.zeros(502)],np.r_[np.array([2]),np.zeros(502)]]
    static_reasons=compact.STATIC_REASONS[1:]
    for chrom in range(1,23):
        rows=[];dose_rows=[]
        for n in range(192+(7 if chrom==1 else 0)):
            snp=f'00_INVENTED_é_χ_chr{chrom}_row{n}';bp=100000+n
            rows.append((chrom,snp,bp));coordinates.append((snp,chrom,bp,'A','C'))
            alleles=('A','C') if n%2==0 else ('C','A')
            db.execute('INSERT INTO source VALUES(?,?,?,?,?,?,?)',(snp,chrom,bp,*alleles,n,3+n*126))
            dose_rows.append(genotype_special[n-192] if n>=192 else base)
        if chrom==6:
            rows.append((chrom,'INVENTED_MHC',25000000))
        if chrom==22:
            for j,reason in enumerate(static_reasons[1:]):
                snp='INVENTED_EXCLUSION_'+reason;bp=200000+j;rows.append((chrom,snp,bp))
                if reason!='MISSING_COORDINATE_MAP_ID':coordinates.append((snp,chrom+(reason=='REFERENCE_COORDINATE_MAP_CONFLICT'),bp,'A' if reason!='INVALID_TARGET_ALLELES' else 'AA','T' if reason=='PALINDROMIC' else 'C'))
                if reason not in ['MISSING_COORDINATE_MAP_ID','MISSING_SOURCE_ID']:
                    al=('G','A') if reason=='ALLELE_CONFLICT_NO_COMPLEMENT_RESCUE' else ('A','C')
                    row0=len(dose_rows);db.execute('INSERT INTO source VALUES(?,?,?,?,?,?,?)',(snp,chrom,bp+(reason=='SOURCE_GRCh37_COORDINATE_CONFLICT'),*al,row0,3+row0*126));dose_rows.append(base)
        ref=refs/f'{chrom}.l2.ldscore.gz'
        with gzip.open(ref,'wt') as f:f.write('CHR SNP BP L2\n');f.writelines(f'{c} {s} {bp} 1\n' for c,s,bp in rows)
        reference_paths.append(str(ref));refhash[str(ref)]=sha(ref);total+=len(rows)
        bed=assets/f'chr{chrom}.bed';bed.write_bytes(b'\x6c\x1b\x01'+b''.join(encode(d) for d in dose_rows));beds[str(chrom)]=str(bed)
    db.commit();db.close()
    with gzip.open(coordinate,'wt') as f:
        f.write('SNP\tCHR\tBP\tA1\tA2\n');f.writelines('\t'.join(map(str,r))+'\n' for r in coordinates)
    refhash[str(coordinate)]=sha(coordinate)
    plan=dict(expected_master_reference_rows=total,reference_body_sha256=refhash,coordinate_map_path=str(coordinate),reference_files_numeric_order=reference_paths)
    r1=assets/'v1';r2=assets/'v2';dirs(r1);dirs(r2);c1.INTENDED_OUTPUT_SHA256.clear();c2.INTENDED_OUTPUT_SHA256.clear()
    static,h=old.prepare_static(plan,r1,source,sha(source),lambda **kw:None)
    sb=new.prepare_static(plan,r2,source,sha(source),lambda **kw:None)
    assert list(compact.static_rows(r2/'qc',sb))==original_rows(static)
    assert list(compact.static_rows(r2/'qc',sb,selected=False))==original_rows(r1/'qc/static_exclusions.jsonl')
    assert sb['logical_sha256']=={'static_selection':sha(static),'static_exclusions':sha(r1/'qc/static_exclusions.jsonl')}
    pass_case('all_static_rows_original_ranks_unicode_leadingzero_alleles_all8_exclusions_exact_v1_dict_and_JSON_SHA')
    j,jh,counts1=old.genotype_qc(r1,static,h,beds,np,lambda **kw:None)
    qb,counts2=new.genotype_qc(r2,sb,beds,np,lambda **kw:None)
    J=list(compact.genotype_rows(r2/'qc',sb,qb));E=list(compact.genotype_rows(r2/'qc',sb,qb,eligible=False))
    assert J==original_rows(j) and E==original_rows(r1/'qc/genotype_exclusions.jsonl') and counts1==counts2
    assert qb['logical_sha256']=={'J_ordered':sha(j),'genotype_exclusions':sha(r1/'qc/genotype_exclusions.jsonl')}
    for a,b in zip(J+E,original_rows(j)+original_rows(r1/'qc/genotype_exclusions.jsonl')):
        for k in ['A1_frequency','MAF']:
            assert a[k] is None and b[k] is None or struct.pack('<d',a[k])==struct.pack('<d',b[k])
    pass_case('all_QC_J_exclusion_rows_exact_integer_statistics_binary64_bits_counts_order_and_v1_JSON_SHA')
    assert qb['QC_rows']==sb['static_rows'] and qb['J_rows']==sum(counts2.values())
    assert total>2048 and all(r['master_rank1']==r['master_rank0']+1 and r['BED_byte_offset']==3+126*r['source_row0'] for r in J)
    pass_case('complete_retention_cross2048_boundaries_original_offsets_and_distinct_J_indices')
    originals={p:sha(p) for p in [r2/'qc/master_static.bin',r2/'qc/master_SNP_utf8.bin',r2/'qc/genotype_QC.bin']}
    def corrupt(role,mutate,binding,consumer):
        path=r2/'qc'/role;before=path.read_bytes();path.write_bytes(mutate(before));changed=dict(binding)
        field={'master_static.bin':'master_sha256','master_SNP_utf8.bin':'SNP_pool_sha256','genotype_QC.bin':'genotype_QC_sha256'}[role];changed[field]=sha(path)
        try:consumer(changed)
        finally:path.write_bytes(before)
    expect('master_new_physical_hash_wrong_logical_SNP_reject',lambda:corrupt('master_SNP_utf8.bin',lambda b:b'X'+b[1:],sb,lambda b:list(compact.static_rows(r2/'qc',b))))
    expect('master_wrong_header_reject',lambda:corrupt('master_static.bin',lambda b:b'X'+b[1:],sb,lambda b:list(compact.static_rows(r2/'qc',b))))
    expect('master_trailing_record_reject',lambda:corrupt('master_static.bin',lambda b:b+b'X',sb,lambda b:list(compact.static_rows(r2/'qc',b))))
    expect('QC_wrong_integer_numerator_even_rehashed_reject',lambda:corrupt('genotype_QC.bin',lambda b:b[:26]+bytes([b[26]^1])+b[27:],qb,lambda b:list(compact.genotype_rows(r2/'qc',sb,b))))
    expect('QC_wrong_original_master_join_reject',lambda:corrupt('genotype_QC.bin',lambda b:b[:16]+struct.pack('<I',999)+b[20:],qb,lambda b:list(compact.genotype_rows(r2/'qc',sb,b))))
    expect('QC_trailing_record_reject',lambda:corrupt('genotype_QC.bin',lambda b:b+b'X',qb,lambda b:list(compact.genotype_rows(r2/'qc',sb,b))))
    assert all(sha(p)==h for p,h in originals.items());pass_case('all_mutation_fixture_bytes_restored_no_original_files_changed')
    writer=compact.IntendedWriter(assets/'tiny_capacity.bin',4)
    try:expect('complete_retention_capacity_failure_before_oversized_write_no_truncation',lambda:writer.write(b'12345'))
    finally:writer.abort()
    assert (assets/'tiny_capacity.bin').stat().st_size==0
    pool_test=assets/'pool_preQC';pool_test.mkdir();pool_writer=compact.StaticWriter(pool_test,1);pool_writer.pool.cap=4
    try:expect('IDpool_whole_stage_fail_before_any_master_or_genotype_row_no_ID_truncation',lambda:pool_writer.append(dict(SNP='INVENTED_TOO_LONG',CHR=1,BP=1,master_rank0=0,master_rank1=1,exclusion='MISSING_SOURCE_ID')))
    finally:pool_writer.abort()
    assert (pool_test/'master_static.bin').stat().st_size==16 and (pool_test/'master_SNP_utf8.bin').stat().st_size==0 and not (pool_test/'genotype_QC.bin').exists()
    sqlpath=assets/'tiny_SQLite_cap.sqlite';capdb=new.configure_db(sqlpath)
    capdb.execute('PRAGMA journal_mode=OFF');capdb.execute('PRAGMA max_page_count=3')
    assert capdb.execute('PRAGMA page_size').fetchone()[0]==4096 and capdb.execute('PRAGMA temp_store').fetchone()[0]==2
    capdb.execute('CREATE TABLE complete(SNP TEXT PRIMARY KEY)')
    try:
        try:
            for i in range(20):capdb.execute('INSERT INTO complete VALUES(?)',(f'INVENTED{i}_'+'x'*1000,))
        except sqlite3.OperationalError as e:assert 'full' in str(e).lower()
        else:raise AssertionError('SQLITE_CAP_EXPECTED_FAIL')
    finally:capdb.close()
    assert sqlpath.stat().st_size<=12288 and not Path(str(sqlpath)+'-journal').exists() and not Path(str(sqlpath)+'-wal').exists()
    pass_case('SQLite_actual_bounded_page_failure_no_disk_journal_temp_no_genotype_QC_credit')
    assert compact.CAPACITY['maximum_planned_QC_bytes']==416822968 and compact.CAPACITY['structural_margin_bytes']==120047944
    assert compact.CAPACITY['maximum_master_and_genotype_fixed_bytes']==64501432
    pass_case('full1290028_structural_bound_includes_DB_pool_exclusions_controls_AppleDouble_metadata_failure_reserve')
    # No curl runs: expose exact actual helper argv with an invented successful
    # transport result; inherited reviewer no-network curlrc witness qualifies -q.
    captured=[]
    def fake_run(command,**kw):
        captured.append(command);Path(command[command.index('--output')+1]).write_bytes(b'INVENTED_ONLY')
        return subprocess.CompletedProcess(command,0,('200\n'+command[-1]+'\n').encode(),b'')
    with patch.object(new.subprocess,'run',fake_run):
        c2.INTENDED_OUTPUT_SHA256.clear();new.curl_file({'curl':'/usr/bin/curl','candidate':{'bytes':288277344}},'https://example.invalid/source',assets/'mock_curl.txt')
    argv=captured[0];assert argv[1]=='-q' and argv.index('-q')==1 and argv[argv.index('--proto')+1]=='=https'
    assert all(x not in argv for x in ['-L','--location','--insecure','-k','--retry','--continue-at','-C'])
    pass_case('actual_curl_helper_first_q_strictHTTPS_no_redirect_resume_retry_invented_transport_only')
    monitor=assets/'monitor.json';mh=save(monitor,dict(command=['WRONG'],status='COMPLETE_REAPED_SOURCE_ONLY',returncode=0,stop_reason=None,cleanup_errors=[],run_error=None,first_teardown_attempt={},process_group_teardown={'remaining_group_members':[]},reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256='P'))
    expect('wrong_frozen_monitor_command_reject_before_worker_result_parse',lambda:controller.consume_worker_result({'output_root':str(assets)},monitor,mh,'P',['EXPECTED']))
    d=assets/'fixed_dirs';dirs(d);f=d/'dummy';f.write_bytes(b'invented');frozen={str(f):sha(f)}
    with patch.object(c2,'freeze_companions',lambda t,m:None):c2.freeze_registered(d,[d,*[p for p in d.rglob('*') if p.is_dir()],f],frozen)
    pass_case('exact_seven_directory_contract_healthy')
    extra=d/'UNEXPECTED_EMPTY';extra.mkdir()
    expect('unexpected_empty_directory_reject_even_dynamic_target_attempt',lambda:c2.freeze_registered(d,[d,*[p for p in d.rglob('*') if p.is_dir()],f],frozen))
    # Scientific functions untouched apart from module namespace imports.
    def funcs(path):
        return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(path).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    assert funcs(P/'scripts/signed_reference_decode_v1.py')==funcs(P/'scripts/signed_reference_decode_v2.py')
    assert funcs(P/'scripts/signed_reference_ldsc_decoderB_v1.py')==funcs(P/'scripts/signed_reference_ldsc_decoderB_v2.py')
    a=funcs(P/'scripts/signed_reference_source_worker_v1.py');b=funcs(P/'scripts/signed_reference_source_worker_v2.py')
    for n in ['line','consume_jsonl','acquire','configure_db','source_format','read_dose']:assert a[n]==b[n],n
    pass_case('unchanged_A_science_and_maintained_B_wrapper_AST_G2_G3_and_source_read_functions')
    # Essential delta integration invokes the actual whole controller with
    # source/locks/runtime/worker/resource explicitly stubbed. Complete compact
    # rows are real invented data; all G2/B receipts remain metadata stubs.
    def whole_controller(name,wrong_command=False,extra_directory=False):
        base=assets/name;base.mkdir();root=base/'canonical_signed_reference_source_v2'
        fixture_plan=dict(output_root=str(root),expected_master_reference_rows=total,source_python=sys.executable,
          source_worker=str(P/'scripts/signed_reference_source_worker_v2.py'),candidate={'bytes':288277344,'upstream_md5':'md5:a7773ab485827b533cb300c76356d76b'},
          resource_contract={'shared_heavy_mutex':str(base/'fake.lock')},raw_transfer_family_mutex=str(base/'fake_raw.lock'),compact_representation=compact.SCHEMA)
        (base/'fake.lock').write_bytes(b'INVENTED');(base/'fake_raw.lock').write_bytes(b'INVENTED')
        pp=base/'plan.json';ph=save(pp,fixture_plan);ap=base/'admission.json';ah=save(ap,{'fixture':'NOT_ACTUAL_ADMISSION'})
        @contextlib.contextmanager
        def fake_lock(*a,**kw):yield 999
        def fake_monitor(plan,command,fds,started,ownership,gates):
            (root/'extracted/members').mkdir();(root/'qc/controls').mkdir()
            c2.INTENDED_OUTPUT_SHA256.clear()
            for target in controller.expected_worker_outputs(root):Path(target).write_bytes(b'INVENTED_METADATA_STUB')
            for f in ['compact_schema.json','master_static.bin','master_SNP_utf8.bin','genotype_QC.bin','reference_index.sqlite','G4_static_seal.json','G4_J_seal.json']:
                shutil.copyfile(r2/'qc'/f,root/'qc'/f)
            shutil.copyfile(source,root/'extracted/source_index.sqlite')
            def replace(path,value):
                Path(path).write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
            ai=dict(bytes=288277344,md5='a7773ab485827b533cb300c76356d76b',sha256=sha(root/'archive/1000G_Phase3_plinkfiles.tgz'))
            replace(root/'qc/G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',candidate=plan['candidate'],archive_identity=ai))
            replace(root/'qc/G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,
              chromosome_variant_counts={str(c):192 for c in range(1,23)},source_index_sha256=sha(root/'extracted/source_index.sqlite')))
            g5records=[]
            for c in range(1,23):
                pref=root/'qc/controls'/f'chr{c}'
                g5records.append(dict(CHR=c,J_count=counts2[str(c)],independent_rows=192,source_window_sha256=sha(str(pref)+'.bed'),
                  B_array_sha256=sha(str(pref)+'.B.npz'),B_receipt_sha256=sha(str(pref)+'.B.receipt.json'),controls=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'],full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
            replace(root/'qc/G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=sha(root/'qc/controls/fixed_G5_selections.json'),chromosomes=g5records,
              maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,finite_reference_503_conditional_only=True,
              population_LD_or_actual_GWAS_mask_adequacy_established=False,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
            mapping={p:sha(p) for p in controller.expected_worker_outputs(root)}
            c2.save(root/'qc/source_worker_result.json',dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=ph,ordinary_output_sha256=mapping,
              archive_identity=ai,source_index_sha256=sha(root/'extracted/source_index.sqlite'),chromosome_J_counts=dict(counts2),J_ordered_sha256=qb['logical_sha256']['J_ordered'],S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
            (root/'worker_stdout.log').write_bytes(b'INVENTED_STUB_NO_WORKER')
            mp=root/'monitor_receipt.json';mh=c2.save(mp,dict(command=['WRONG'] if wrong_command else command,status='COMPLETE_REAPED_SOURCE_ONLY',returncode=0,stop_reason=None,
              cleanup_errors=[],run_error=None,first_teardown_attempt={},process_group_teardown={'remaining_group_members':[]},reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256=ph))
            if extra_directory:(root/'UNEXPECTED_EMPTY').mkdir()
            return mp,mh
        with patch.object(controller,'exclusive_heavy_lock',fake_lock),patch.object(controller,'shared_raw_family_lock',fake_lock),\
             patch.object(controller,'identity_gate',lambda *a:None),patch.object(controller,'runtime_profile_gate',lambda *a:None),\
             patch.object(controller,'resource_snapshot',lambda *a,**kw:{}),patch.object(controller,'monitor_worker',fake_monitor),\
             patch.object(controller,'cleanup',lambda *a:None),patch.object(c2,'resource_snapshot',lambda *a,**kw:{}):
            if wrong_command or extra_directory:
                expect(name,lambda:controller.execute(pp,ph,ap,ah))
                assert (root/'PENDING.json').exists() and not (root/'terminal_seal.json').exists()
            else:
                controller.execute(pp,ph,ap,ah)
                assert not (root/'PENDING.json').exists() and (root/'terminal_seal.json').exists()
                pass_case(name)
    whole_controller('whole_source_controller_stub_compact_full_reconstruction_and_private_Terminal2_healthy')
    whole_controller('whole_source_controller_stub_wrong_frozen_monitor_command_reject_PENDING',wrong_command=True)
    whole_controller('whole_source_controller_stub_unexpected_empty_directory_reject_PENDING',extra_directory=True)
    meter_root=assets/'whole_source_controller_stub_compact_full_reconstruction_and_private_Terminal2_healthy/canonical_signed_reference_source_v2'
    contract=dict(source_epoch_namespaces=[str(meter_root.with_name('canonical_signed_reference_source_v1')),str(meter_root)],
      source_aggregate_cap_bytes=3<<30,archive_and_all_failed_partial_cap_bytes=512<<20,extracted_all_partial_and_failed_cap_bytes=2<<30,
      QC_controls_metadata_and_all_failed_output_cap_bytes=512<<20,global_ceiling_bytes=300<<30,closed_sensitivity_cap_bytes=256<<20,
      internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,SSD_launch_free_required_bytes=8<<30,whole_stage_monotonic_deadline_seconds=7200)
    mp=dict(output_root=str(meter_root),resource_contract=contract,compact_representation=compact.SCHEMA,expected_master_reference_rows=total,closed_sensitivity_namespace=str(assets/'CLOSED_INVENTED'))
    real_tree=c2.full_tree_bytes
    def old_failed_epoch_meter(path):
        p=Path(path)
        if p==meter_root.with_name('canonical_signed_reference_source_v1'):return 600<<20
        return real_tree(p)
    with patch.object(c2.shutil,'disk_usage',lambda p:SimpleNamespace(free=20<<30)),patch.object(c2,'full_tree_bytes',old_failed_epoch_meter):
        expect('prior_epoch_failed_QC_bytes_count_against_unchanged_family512MiB_no_quota_reuse',lambda:c2.resource_snapshot(mp,c2.time.monotonic()))
    receipt=dict(status='AUTHOR_FOCUSED_COMPACT_EQUIVALENCE_AND_OPERATIONAL_DELTA_PASS_NO_EXECUTION_ADMISSION',cases=cases,case_count=len(cases),
      runtime_python=sys.version,numpy=np.__version__,invented_master_rows=total,static_binding=sb,QC_binding=qb,
      exact_four_original_logical_streams=True,exact_binary64_frequency_and_MAF_bits=True,structural_capacity=compact.CAPACITY,
      no_source_actual_capacity_success_claim=True,actual_source_body_network_production_workers_or_mutex_operations=0,
      S_simulations_run=0,C_calibration_credit=0,estimator_calls=0)
    h=save(OUT/'delta_control_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],case_count=len(cases),receipt_sha256=h)))


if __name__=='__main__':main()
