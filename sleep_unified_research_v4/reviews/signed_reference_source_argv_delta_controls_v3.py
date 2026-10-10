"""Two metadata-stub argv-fence controls; no compact/source/runtime replay."""
import ast
import contextlib
import hashlib
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';sys.path.insert(0,str(S))
import signed_reference_source_controller_v3 as controller
import signed_reference_source_common_v3 as common
from signed_reference_compact_v3 import SCHEMA,SCHEMA_SHA256


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,value):
    b=(json.dumps(value,sort_keys=True,indent=2)+'\n').encode();h=hashlib.sha256(b).hexdigest()
    with Path(p).open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
    assert sha(p)==h;return h
def funcs(p):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(p).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}


def main():
    assert sys.version_info[:3]==(3,11,11)
    out=Path(sys.argv[1]);out.mkdir();results=[]
    source_function_equivalence={}
    for stem in ['signed_reference_source_worker','signed_reference_compact','signed_reference_ldsc_decoderB']:
        assert funcs(S/(stem+'_v2.py'))==funcs(S/(stem+'_v3.py')),stem
        source_function_equivalence[stem]='IDENTICAL_FUNCTION_AST_IMPORT_NAMESPACE_ONLY'
    a=funcs(S/'signed_reference_source_controller_v2.py');b=funcs(S/'signed_reference_source_controller_v3.py')
    for name in a:
        if name!='execute':assert a[name]==b[name],name
    a=funcs(S/'signed_reference_source_common_v2.py');b=funcs(S/'signed_reference_source_common_v3.py')
    for name in a:
        if name!='resource_snapshot':assert a[name]==b[name],name
    capacity=SCHEMA['capacity'];assert capacity['maximum_planned_QC_bytes']==416822968 and SCHEMA_SHA256=='26b4fdc33d6ecf81a0c26a18eb59249da1ff0c86587110f61f00a86be28ac074'
    # Prove the only executable controller AST change is moving original tuple
    # capture before the callback and passing a separate list to that callback.
    text=(S/'signed_reference_source_controller_v2.py').read_text()
    old="                    monitor,msha=monitor_worker(plan,command,(heavy_fd,raw_fd),started,ownership,dict(identity=identity,plan_sha=plan_sha))\n                    command=tuple(command)"
    new="                    command=tuple(command)\n                    monitor,msha=monitor_worker(plan,list(command),(heavy_fd,raw_fd),started,ownership,dict(identity=identity,plan_sha=plan_sha))"
    assert old in text
    intended=ast.parse(text.replace(old,new))
    desired=next(n for n in intended.body if isinstance(n,ast.FunctionDef) and n.name=='execute')
    assert ast.dump(desired,include_attributes=False)==funcs(S/'signed_reference_source_controller_v3.py')['execute']
    @contextlib.contextmanager
    def fake_lock(*a,**kw):yield 999
    for mutate in [False,True]:
        name='MUTABLE_ARGV_NEGATIVE' if mutate else 'HEALTHY';base=out/name;base.mkdir();root=base/'canonical_signed_reference_source_v3'
        lock=base/'INVENTED.lock';lock.write_bytes(b'INVENTED_NO_FLOCK')
        plan=dict(output_root=str(root),expected_master_reference_rows=4224,source_python='/INVENTED_READ_ONLY_PYTHON',
          source_worker=str(S/'signed_reference_source_worker_v3.py'),candidate=dict(bytes=288277344,upstream_md5='md5:a7773ab485827b533cb300c76356d76b'),
          resource_contract={'shared_heavy_mutex':str(lock)},raw_transfer_family_mutex=str(lock),compact_representation=SCHEMA)
        pp=base/'STUB_PLAN_NOT_ACTUAL.json';ph=save(pp,plan);ap=base/'STUB_ADMISSION_NOT_ACTUAL.json';ah=save(ap,{'not_actual_admission':True})
        original=[];seen_expected=[];seen_callback=[]
        real_consume=controller.consume_worker_result
        def traced_consume(*args,**kwargs):
            expected=args[4];assert isinstance(expected,tuple);seen_expected.append(expected)
            assert expected==original[0],'ORIGINAL_EXPECTED_TUPLE_REBOUND'
            return real_consume(*args,**kwargs)
        def fake_monitor(plan,command,fds,started,ownership,gates):
            assert isinstance(command,list);original.append(tuple(command));seen_callback.append(command)
            if mutate:command[0]='INVENTED_MUTATED_EXECUTABLE'
            assert original[0][0]=='/INVENTED_READ_ONLY_PYTHON'
            (root/'extracted/members').mkdir();(root/'qc/controls').mkdir();common.INTENDED_OUTPUT_SHA256.clear()
            for path in controller.expected_worker_outputs(root):Path(path).write_bytes(b'INVENTED_NO_SOURCE_OR_COMPACT_STREAM')
            def replace(p,v):Path(p).write_text(json.dumps(v,sort_keys=True,indent=2)+'\n')
            replace(root/'qc/compact_schema.json',SCHEMA)
            sb=dict(schema_sha256=SCHEMA_SHA256,master_rows=4224,static_rows=4224,master_sha256=sha(root/'qc/master_static.bin'),
              SNP_pool_sha256=sha(root/'qc/master_SNP_utf8.bin'),logical_sha256={'static_selection':'1'*64,'static_exclusions':'2'*64})
            qb=dict(schema_sha256=SCHEMA_SHA256,QC_rows=4224,J_rows=4224,genotype_QC_sha256=sha(root/'qc/genotype_QC.bin'),
              logical_sha256={'J_ordered':'3'*64,'genotype_exclusions':'4'*64})
            replace(root/'qc/G4_static_seal.json',dict(status='FROZEN_BEFORE_GENOTYPE_QC',master_reference_rows=4224,static_rows=4224,
              compact_schema_file_sha256=sha(root/'qc/compact_schema.json'),compact_static=sb,exclusion_counts={},static_ordered_sha256='1'*64,static_exclusions_sha256='2'*64))
            counts={str(c):192 for c in range(1,23)}
            replace(root/'qc/G4_J_seal.json',dict(status='REFERENCE_ONLY_J_FROZEN',compact_genotype=qb,J_count=4224,J_ordered_sha256='3'*64,
              chromosome_J_counts=counts,genotype_exclusion_counts={},genotype_exclusions_sha256='4'*64,no_MAF_threshold=True,no_GWAS_outcomes=True,no_native_mask_selection=True))
            ai=dict(bytes=288277344,md5='a7773ab485827b533cb300c76356d76b',sha256=sha(root/'archive/1000G_Phase3_plinkfiles.tgz'))
            replace(root/'qc/G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',candidate=plan['candidate'],archive_identity=ai))
            replace(root/'qc/G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,
              chromosome_variant_counts=counts,source_index_sha256=sha(root/'extracted/source_index.sqlite')))
            records=[]
            for c in range(1,23):
                pref=root/'qc/controls'/f'chr{c}'
                records.append(dict(CHR=c,J_count=192,independent_rows=192,source_window_sha256=sha(str(pref)+'.bed'),B_array_sha256=sha(str(pref)+'.B.npz'),
                  B_receipt_sha256=sha(str(pref)+'.B.receipt.json'),controls=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'],
                  full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
            replace(root/'qc/G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=sha(root/'qc/controls/fixed_G5_selections.json'),chromosomes=records,
              maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,finite_reference_503_conditional_only=True,population_LD_or_actual_GWAS_mask_adequacy_established=False,
              S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
            mapping={p:sha(p) for p in controller.expected_worker_outputs(root)}
            common.save(root/'qc/source_worker_result.json',dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=ph,ordinary_output_sha256=mapping,
              archive_identity=ai,source_index_sha256=sha(root/'extracted/source_index.sqlite'),chromosome_J_counts=counts,J_ordered_sha256='3'*64,
              S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
            (root/'worker_stdout.log').write_bytes(b'INVENTED_NO_WORKER');mp=root/'monitor_receipt.json'
            mh=common.save(mp,dict(command=command,status='COMPLETE_REAPED_SOURCE_ONLY',returncode=0,stop_reason=None,cleanup_errors=[],run_error=None,
              first_teardown_attempt={},process_group_teardown={'remaining_group_members':[]},reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256=ph))
            return mp,mh
        def static_stub(*a,**kw):return iter([])
        def genotype_stub(*a,**kw):return iter([] if kw.get('eligible') is False else [{'CHR':c} for c in range(1,23) for _ in range(192)])
        with patch.object(controller,'exclusive_heavy_lock',fake_lock),patch.object(controller,'shared_raw_family_lock',fake_lock),\
             patch.object(controller,'identity_gate',lambda *a:None),patch.object(controller,'runtime_profile_gate',lambda *a:None),\
             patch.object(controller,'resource_snapshot',lambda *a,**kw:{}),patch.object(common,'resource_snapshot',lambda *a,**kw:{}),\
             patch.object(controller,'monitor_worker',fake_monitor),patch.object(controller,'cleanup',lambda *a:None),\
             patch.object(controller,'consume_worker_result',traced_consume),patch.object(controller,'static_rows',static_stub),patch.object(controller,'genotype_rows',genotype_stub):
            if mutate:
                try:controller.execute(pp,ph,ap,ah)
                except RuntimeError as e:assert 'SOURCE_EXACT_CONTROLLER_DERIVED_COMMAND_REQUIRED' in str(e)
                else:raise AssertionError('MUTABLE_ARGV_MUST_REJECT')
                assert (root/'PENDING.json').exists() and not (root/'terminal_seal.json').exists() and (root/'source_execution_receipt.json.failure.json').exists()
                assert len(seen_expected)==1 and seen_callback[0][0]=='INVENTED_MUTATED_EXECUTABLE'
            else:
                controller.execute(pp,ph,ap,ah)
                assert not (root/'PENDING.json').exists() and (root/'terminal_seal.json').exists() and len(seen_expected)>=2
                assert all(x is seen_expected[0] for x in seen_expected),'ORIGINAL_TUPLE_OBJECT_NOT_REUSED'
        results.append(dict(case=name,status='PASS',initial_final_original_expected_tuple_unchanged=True,
          callback_received_separate_mutable_list=True,actual_compact_streams_and_all_source_runtime_locks_workers_resources='STUBBED_NO_CREDIT'))
    h=save(out/'argv_delta_receipt.json',dict(status='AUTHOR_TWO_FOCUSED_ARGV_FENCE_CONTROLS_PASS_NO_EXECUTION_ADMISSION',case_count=2,cases=results,
      unchanged_import_only_scientific_function_AST=source_function_equivalence,unchanged_compact_schema_sha256=SCHEMA_SHA256,
      only_controller_execute_AST_delta='FREEZE_BEFORE_CALLBACK_AND_PASS_LISTCOPY',
      source_epoch_namespaces_v1_v2_v3_in_original_family_budget=True,old23_42_suites_compact_math_streams_and_runtime_census_repeated=0,
      actual_source_archive_genotype_GWAS_reference_body_operations=0,actual_transfers_production_workers_or_mutex_operations=0,
      S_simulations_run=0,C_calibration_credit=0,estimator_calls=0))
    print(json.dumps(dict(status='AUTHOR_ARGV_DELTA_PASS',case_count=2,sha256=h)))


if __name__=='__main__':main()
