"""Prepare only: preserve7 successes and snoring; replay remaining27 once."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

from extension_replay_common_v3 import sha, write_new
from core_large35_checkpoint_adoption_v1 import regular, qc, semantic_snoring, adoption_gate, FIRST8

P = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OLD = SSD/'core_pipeline/large35_bounded_replay_v3'
OUT = SSD/'core_pipeline/large35_checkpoint_continuation_v1'
OLD_PLAN_SHA = '3e8c3c7bca26b00ed337c816a723d8a42475473667d4a9c969348045fcc9a2f4'
METER_SHA = '6fec56a75ad71d841379436a461ac63d0d8cbcbe0ac62cf4aa960ead6a783d88'


def main():
    old_path = regular(OLD/'core_large35_bounded_replay_plan_v3.json')
    if sha(old_path) != OLD_PLAN_SHA:raise RuntimeError('EXACT_PRESERVED_FAILED_CORE_PLAN_REQUIRED')
    old = json.loads(old_path.read_text())
    execution_path = regular(OLD/'core_large35_execution_receipt_v3.json')
    execution_sha = sha(execution_path); failed = json.loads(execution_path.read_text())
    if (failed['status'] != 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' or failed['owned_cleanup_verified'] is not True
        or failed['termination_requests'] or [m['trait'] for m in failed['completed_members']] != FIRST8[:7]
        or len(failed['worker_receipt_sha256']) != 38):raise RuntimeError('EXACT_ACTUAL7_SUCCESS_SNORING_STOP_REQUIRED')
    if OUT.exists() or OUT.is_symlink():raise RuntimeError('NEW_CONTINUATION_NAMESPACE_MUST_BE_UNUSED')
    def relocate(value):
        if isinstance(value, dict):return {k: relocate(v) for k,v in value.items()}
        if isinstance(value, list):return [relocate(v) for v in value]
        if isinstance(value, str) and value.startswith(str(OLD)+'/'):return str(OUT)+value[len(str(OLD)):]
        return value
    members = [relocate(m) for m in old['members'][8:]]
    metadata = {str(old_path):OLD_PLAN_SHA,str(execution_path):execution_sha}
    def bind(path, expected=None):
        value=sha(regular(path))
        if expected is not None and value!=expected:raise RuntimeError('ACTUAL_SMALL_CHECKPOINT_IDENTITY_CHANGED: '+str(path))
        metadata[str(path)]=value
        return value
    bind(OLD/'core_large35_pending_v3.json')
    adopted=[]
    for member, result in zip(old['members'][:7],failed['completed_members']):
        adopted.append(dict(trait_id=member['trait_id'],origin='PREDECESSOR_COMPLETED_EXACT_CONTENT_QC',member=member,output_sha256=result['output_sha256']))
        for path,digest in result['output_sha256'].items():
            if not path.endswith(('.harmonized.tsv.gz','.sumstats.gz')):bind(path,digest)
    for path,digest in failed['worker_receipt_sha256'].items():bind(path,digest)
    member=old['members'][7];spool=Path(member['ephemeral_spool'])
    candidate=spool/'global_harmonization_candidate_receipt.json';columns=spool/'columns/column_preparation_manifest.json'
    for path,digest in failed['generated_metadata_sha256']['snoring'].items():bind(path,digest)
    for key in ['source_gate_receipt','harmonization_qc','comparison_receipt']:bind(member[key])
    original_qc=member['original_design']['historical_QC'];bind(original_qc['path'],original_qc['sha256'])
    failed_worker=OLD/'receipts_v4/snoring__compare.worker.json';bind(failed_worker)
    worker=json.loads(failed_worker.read_text());comparison=json.loads(Path(member['comparison_receipt']).read_text())
    if (worker['status']!='WORKER_FAILED_PRESERVED' or worker['returncode']!=1 or worker['stop_reason'] is not None
        or worker['metadata_errors'] or worker['owned_cleanup_verified'] is not True
        or worker['process_group_teardown']['remaining_group_members']
        or worker['output_sha256'].get(member['comparison_receipt'])!=metadata[member['comparison_receipt']]
        or comparison['plan_sha256']!=OLD_PLAN_SHA):raise RuntimeError('EXACT_REAPED_LABEL_COMPARISON_FAILURE_REQUIRED')
    observed, steps=qc(member['harmonization_qc']);column_data=json.loads(columns.read_text())
    semantic=semantic_snoring('snoring',original_qc,observed,steps,column_data,comparison)
    panel=P.parent/'config/analysis_panel.tsv';bind(panel)
    with panel.open() as f:row=next(r for r in csv.DictReader(f,delimiter='\t') if r['trait_id']=='snoring')
    if row['type']!='binary' or float(row['ncase'])!=152302 or float(row['ncontrol'])!=256015:
        raise RuntimeError('ACTUAL_CONSTANT_EFFECTIVE_N_CONFIG_CHANGED')
    current=P.parent/'scripts/01_harmonize.py';bind(current,old['dependencies_sha256'][str(current)])
    historical=subprocess.run(['git','show','546a0464^:scripts/01_harmonize.py'],cwd=P.parent,capture_output=True,check=True).stdout
    evidence=dict(historical_revision='546a0464^',historical_code_sha256=hashlib.sha256(historical).hexdigest(),
                  current_code_path=str(current),current_code_sha256=sha(current),
                  predicate='N.isna() | (N <= 0) | (N < 0.5 * configured_effective_N)',
                  historical_branch='total_effective_n = effective_n(metadata_ncase, metadata_ncontrol); out[N] = total_effective_n',
                  current_branch='n_reference = effective_n(metadata_ncase, metadata_ncontrol); out[N] = n_reference',
                  configured_counts={'ncase':152302,'ncontrol':256015,'n_total':408317},
                  historical_per_trait_execution_attestation_claimed=False,
                  other_per_SNP_N_branches_not_adjudicated=True)
    proof_path=P/'statistical_validation/snoring_constant_N_QC_semantic_adjudication_v1.json'
    proof=dict(schema='actual_snoring_constant_N_label_adjudication_v1',semantic_decision=semantic,
               evidence_sha256=dict(metadata),predicate_evidence=evidence,
               status='ACTUAL_MATCHED_CONTENT_WITH_SCOPED_SEMANTIC_QC_LABEL_QUALIFICATION',
               source_body_reads=0,output_body_reads=0,workers_launched=0,rerun_matched_content=False)
    write_new(proof_path,proof);bind(proof_path)
    outputs={member['harmonized']:comparison['harmonized_comparison']['compressed_sha256'][0],
             member['munged']:comparison['munged_comparison']['compressed_sha256_new']}
    outputs.update({path:metadata[path] for path in [member['source_gate_receipt'],member['harmonization_qc'],member['comparison_receipt'],str(candidate),str(columns)]})
    adopted.append(dict(trait_id='snoring',origin='SNORING_MATCHED_CONTENT_SCOPED_SEMANTIC_QC_ADJUDICATION',member=member,
                        output_sha256=outputs,candidate_path=str(candidate),columns_path=str(columns),
                        failed_spool_preserved=True))
    plan=dict(old)
    dependencies=dict(old['dependencies_sha256'])
    for name in ['core_large35_checkpoint_adoption_v1.py','core_large35_checkpoint_continuation_v1.py','prepare_core_large35_checkpoint_continuation_v1.py','108_acquire_validation_raw_sources_v4.py']:
        file=P/'scripts'/name;dependencies[str(file)]=sha(regular(file))
    if dependencies[str(P/'scripts/108_acquire_validation_raw_sources_v4.py')]!=METER_SHA:raise RuntimeError('EXACT_REVIEWED_NOFOLLOW108_METER_REQUIRED')
    dependencies.update(metadata)
    plan.update(schema='actual_core_large35_checkpoint_continuation_plan_v1',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                predecessor_plan=str(old_path),predecessor_plan_sha256=OLD_PLAN_SHA,predecessor_execution_receipt=str(execution_path),
                predecessor_failed_family_reclassified=False,private_namespace=str(OUT),members=members,member_count=27,
                adopted_members=adopted,adopted_worker_receipt_sha256=failed['worker_receipt_sha256'],
                original35_order=[m['trait_id'] for m in old['members']],checkpoint_metadata_sha256=metadata,
                snoring_adjudication=str(proof_path),dependencies_sha256=dependencies,
                deadline_anchor_utc=old['prepared_utc'],
                deadline_qualification='Conservative original preparation UTC precedes actual original start; same96h budget, never reset.',
                meter_source=str(P/'scripts/108_acquire_validation_raw_sources_v4.py'),meter_source_sha256=METER_SHA,
                cumulative_capacity='Original16GiB across all core_pipeline v1/v2/v3/continuation outputs and failed spools; global300GiB unchanged.',
                preparation_identity_qualification='Only small checkpoint/code/config metadata freshly hashed. Inherited body/output/reference/runtime identities not freshly verified here; full current dependency/runtime/retained-output checks remain mandatory at execution.',
                expected_new_worker_count=135,adopted_result_count=8,combined_completion_count=35,
                source_body_reads_in_preparation=0,output_body_reads_in_preparation=0,workers_launched=0,execution_admitted=False)
    if adoption_gate(plan,current_content_hashes=False)['adopted_members']!=8:raise RuntimeError('ADOPTION_PREPARATION_FAILED')
    for m in members:
        for folder in ['prefilter','harmonized','munged','receipts']:(OUT/m['trait_id']/folder).mkdir(parents=True,exist_ok=False)
    for folder in ['receipts_v4','logs_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=False)
    path=OUT/'core_large35_checkpoint_continuation_plan_v1.json';write_new(path,plan)
    print(json.dumps(dict(plan=str(path),plan_sha256=sha(path),snoring_adjudication=str(proof_path),snoring_adjudication_sha256=sha(proof_path),adopted=8,remaining=27,new_workers=135,workers_launched=0,body_reads=0),indent=2))


if __name__=='__main__':main()
