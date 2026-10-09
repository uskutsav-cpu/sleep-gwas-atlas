#!/usr/bin/env python3
"""Consolidate technical evidence; no manuscript or association narrative."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

PACKAGE=Path(__file__).resolve().parents[1];ROOT=PACKAGE.parent

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def save(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2);f.write('\n')

def textfile(path,value):
    with path.open('x') as f:f.write(value)

def table(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def load(relative):return json.loads((PACKAGE/relative).read_text())

def main():
    for trait in ['breast_cancer','ovarian_cancer']:
        if load('logs/'+trait+'_acquisition_receipt_v3.json')['status']!='EXACT_HISTORICAL_SOURCE_REACQUIRED':raise RuntimeError('TWO_EXACT_CORE_CONTAINERS_REQUIRED')
        if load('logs/'+trait+'_materialization_receipt_v3.json')['status']!='NATIVE_RAW_MATERIALIZATION_EXACT_CONTENT':raise RuntimeError('TWO_RAW_CONTENT_REPRODUCTIONS_REQUIRED')
    if load('reviews/mvp_preprocessing_stream_review_v3.json')['status']!='INDEPENDENT_STREAMING_PREPROCESSING_PASS':raise RuntimeError('INDEPENDENT_MVP_REVIEW_REQUIRED')
    seals=[]
    for name,manifest in [('sleep_unified_research_v1','hashes/FILES_SHA256.tsv'),('sleep_unified_research_v2','hashes/ARTIFACT_SHA256_v2.tsv')]:
        for r in read(ROOT/name/manifest):
            p=ROOT/name/r['path']
            if sha(p)!=r['sha256'] or p.stat().st_size!=int(r['bytes']):raise RuntimeError('EARLIER_SEAL_CHANGED '+str(p))
            seals.append({'package':name,'path':r['path'],'sha256':r['sha256'],'status':'UNCHANGED'})
    if len(seals)!=278:raise RuntimeError('EARLIER_SEAL_CARDINALITY_CHANGED')
    table(PACKAGE/'tables/EARLIER_SEAL_VERIFICATION_V3.tsv',seals)
    v2=ROOT/'sleep_unified_research_v2'
    current=read(v2/'tables/CURRENT_CORE_SOURCE_STATUS.tsv'); updates=read(v2/'tables/SOURCE_RECOVERY_UPDATES.tsv')
    acquired={}
    for trait in ['breast_cancer','prostate_cancer','ovarian_cancer']:
        acquired[trait]=load('logs/'+trait+'_acquisition_receipt_v3.json')
    for r in current:
        r.update(public_candidate_path='',public_candidate_sha256='',public_candidate_bytes='',public_candidate_status='')
        if r['kind']=='core_source_archive' and r['trait_id'] in acquired:
            a=acquired[r['trait_id']];r.update(public_candidate_path=a['path'],public_candidate_sha256=a['actual_sha256'],public_candidate_bytes=a['actual_bytes'],public_candidate_status=a['status'])
            if a['status']=='EXACT_HISTORICAL_SOURCE_REACQUIRED' and a['actual_sha256']==r['expected_sha256']:
                r.update(current_path=a['path'],current_bytes=a['actual_bytes'],current_sha256=a['actual_sha256'],current_status='EXACT_HISTORICAL_HASH_AVAILABLE',proof='logs/'+r['trait_id']+'_acquisition_receipt_v3.json')
            else:r.update(current_status='UNRECOVERED_EXACT_CONTAINER_PUBLIC_HASH_MISMATCH',proof='logs/'+r['trait_id']+'_acquisition_receipt_v3.json')
    table(PACKAGE/'tables/CURRENT_CORE_SOURCE_STATUS_V3.tsv',current)
    for r in updates:
        if r['trait_id'] in acquired and acquired[r['trait_id']]['status']=='EXACT_HISTORICAL_SOURCE_REACQUIRED':
            a=acquired[r['trait_id']]
            if a['actual_sha256']!=r['expected_sha256']:raise RuntimeError('INITIAL_MISSING_LEDGER_HASH_MISMATCH')
            r.update(resolved_path=a['path'],resolved_bytes=a['actual_bytes'],actual_sha256=a['actual_sha256'],recovery_method='EXACT_HISTORICAL_PUBLIC_FILE_REACQUIRED_V3',current_status='EXACT_EXPECTED_SOURCE_BYTES_AVAILABLE_SEPARATE_LOCATION')
    table(PACKAGE/'tables/CUMULATIVE_MISSING_SOURCE_RECOVERY_V3.tsv',updates)
    native=[]
    for trait in ['breast_cancer','ovarian_cancer']:
        d=load('logs/'+trait+'_materialization_receipt_v3.json'); old=d['historical_decompressed'];new=d['native_decompressed']
        native.append({'trait_id':trait,'status':d['status'],'old_decompressed_sha256':old['decompressed_sha256'],'new_decompressed_sha256':new['decompressed_sha256'],
                       'data_rows':new['newline_count_including_header']-1,'decompressed_bytes':new['decompressed_bytes'],
                       'gzip_bytes_identical':d['gzip_bytes_identical'],'numeric_or_content_difference':'NONE;EXACT_DECOMPRESSED_BYTES',
                       'serialization_difference':'GZIP_SERIALIZATION_DIFF;INDEPENDENT_HEADER_DEFLATE_REVIEW_REQUIRED' if not d['gzip_bytes_identical'] else 'NONE',
                       'input_hashes_unchanged':d['input_hashes_unchanged'],'peak_observed_worker_RSS_bytes':d['worker_peak_observed_RSS_bytes'],
                       'elapsed_seconds':d['elapsed_seconds_including_comparison'],'full_GWAS_QC_or_estimator_reproduced':False,'proof':'logs/'+trait+'_materialization_receipt_v3.json'})
    table(PACKAGE/'tables/NATIVE_CANCER_RAW_CONTENT_COMPARISON_V3.tsv',native)
    rejected=load('logs/rejected_prostate_metadata_diagnostic_v3.json')
    table(PACKAGE/'tables/REJECTED_PROSTATE_SOURCE_COUNT_COMPARISON_V3.tsv',[{'metric':k,**v,'source_still_rejected':True,'cause_resolved':False} for k,v in rejected['comparison_to_historical'].items()])
    head=load('logs/extension_version_availability_v3.json')
    table(PACKAGE/'tables/EXTENSION_LOCKED_VERSION_AVAILABILITY_V3.tsv',[{'trait_id':r['extension_trait_id'],'expected_bytes':r['expected_bytes'],
         'expected_checksum':r['expected_checksum'],'expected_version_id':r['expected_version_id'],'status':r['status'],
         'current_body_checksum_verified':False,'header_proof':r['header_path'],'header_sha256':r['header_sha256']} for r in head['results']])
    full={'recorded_utc':datetime.now(timezone.utc).isoformat(),'internal_free_bytes':shutil.disk_usage(ROOT).free,
          'ssd_free_bytes':shutil.disk_usage('/Volumes/Extreme SSD').free,'internal_minimum_bytes':3*1024**3,'ssd_minimum_bytes':5*1024**3,
          'memory':subprocess.check_output(['sysctl','hw.memsize','vm.swapusage'],text=True),
          'native_runner_sha256':sha(ROOT/'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'),
          'full_native_job_launched':False,'original_sealed_native_runner_or_guard_changed':False}
    full['resource_gate_pass']=full['internal_free_bytes']>=full['internal_minimum_bytes'] and full['ssd_free_bytes']>=full['ssd_minimum_bytes']
    save(PACKAGE/'logs/native_resource_preflight_v3.json',full)
    prior=[ROOT/'sleep_unified_research_v1/logs/native_resource_preflight_20261009T045116Z.json',ROOT/'sleep_unified_research_v2/logs/native_resource_preflight_v2.json',PACKAGE/'logs/native_resource_preflight_v3.json']
    audit={'same_full_native_resource_condition_in_three_consecutive_goal_turns':all(not json.loads(p.read_text())['resource_gate_pass'] for p in prior),
           'proofs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'preflight':json.loads(p.read_text())} for p in prior],
           'useful_completed_work_this_turn':['two exact original container recoveries','two unchanged native raw projections','MVP whole-source preprocessing and independent numerical verification','100 locked source HEAD identities','rejected prostate forensic metadata'],
           'remaining_scientific_jobs_require_unchanged_full_native_guard':True,
           'massive_extension_body_mirror_not_started':'Native preprocessing/estimator resource gate fails; exact immutable remote versions remain available. A 212GiB mirror alone would not complete a scientific analysis.',
           'no_primary_confirmatory_local_or_molecular_experiment_admitted':True,'required_external_change':'Free at least3GiB internally or provide a suitable compute host; resolve provenance/CI/covariance human questions before inferential admission',
           'goal_not_complete':True}
    save(PACKAGE/'manifests/full_native_blocking_audit_v3.json',audit)
    mvp=load('logs/mvp_preprocessing_worker_receipt_v3.json');mvp_review=load('reviews/mvp_preprocessing_stream_review_v3.json')
    state={'recorded_utc':datetime.now(timezone.utc).isoformat(),'initial_missing_ledger_entries':len(updates),
           'cumulative_missing_entries_resolved':sum(r['current_status']!='UNRECOVERED' for r in updates),
           'cumulative_unique_missing_hashes_recovered':len({r['actual_sha256'] for r in updates if r['actual_sha256']}),
           'raw_source_hashes_available':sum(r['kind']=='core_raw' and r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in current),
           'exact_core_container_hashes_available':sum(r['kind']=='core_source_archive' and r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in current),
           'remaining_exact_container_unrecovered':['prostate_cancer'],'native_raw_projections_reproduced_this_version':2,
           'native_lipid_prefilters_reproduced_previous_version':3,'native_core_rg_controls_total':1,
           'native_MVP_preprocessing_source_rows':mvp['counts']['source_rows'],'native_MVP_retained_HM3_rows':mvp['counts']['output_rows'],
           'independent_MVP_checks':mvp_review['check_count'],'extension_locked_versions_HEAD_matched':head['matched_HEAD_versions'],
           'extension_body_bytes_downloaded_this_version':0,'complete_core_396_native_reproduction':False,
           'complete_extension_1200_native_reproduction':False,'complete_validation_41_native_reproduction':False,
           'both_trait_independent_successes':0,'new_pair_outcomes_estimated':0,'new_primary_question_admitted':False,
           'shared_covariance_calibrated':False,'local_or_molecular_experiment_admitted':False,
           'HUMAN_MANUSCRIPT_AUTHORING_READY':False,'prior_sealed_files_verified_unchanged':len(seals),
           'full_native_resource_gate_pass':full['resource_gate_pass'],'scientific_goal_complete':False,
           'source_status_scope':'v1/v2 verified source hashes carried forward; newly recovered containers and projected raw inputs rehashed before/after current execution'}
    save(PACKAGE/'tables/CURRENT_EVIDENCE_STATUS_V3.json',state)
    gates=[('Public clinical definition','PASS_QUALIFIED','GIA/CIPHER clinical insomnia; not exact UKB symptom equivalence'),
           ('Source identity','PASS','575504178bytes and acquisition MD5/SHA confirmed independently'),
           ('Coordinate and allele concordance','PASS_PREPROCESSING_ONLY','826027retained; independently verified reference coordinates/alleles'),
           ('Upstream uncertainty and effect semantics','UNRESOLVED','95% log-scale Wald CI and conversion/effect semantics require source-specific confirmation'),
           ('Source h2 and intercept','NOT_ESTIMATED_RESOURCE_GATED','Unchanged full-native3GiB internal guard'),
           ('Power','NOT_PASSED','No source h2 or prospective pair uncertainty yet'),
           ('Participant independence','QUALIFIED_UNKNOWN_INTERSECTION','Named MVP/UKB/FinnGen cohorts; individual cross-enrollment unknown'),
           ('Pair admission','ZERO','Ten metadata-selected transport candidates remain conditional'),
           ('Confirmatory novelty','NO_GO_CURRENT','No new primary question passes the locked feasibility matrix')]
    table(PACKAGE/'tables/MVP_DOWNSTREAM_GATES_V3.tsv',[{'gate':g,'status':s,'qualification':q} for g,s,q in gates])
    table(PACKAGE/'tables/PHASE_COMPLETION_V3.tsv',[{'phase':p,'status':s,'evidence_summary':e} for p,s,e in [
        (0,'PARTIAL_SOURCE_PROVENANCE','45/45 raw;42/43 exactcontainers;publicprostate rejected withrealrow differences'),
        (1,'PARTIAL','1nativecorecontrol;3lipidprefilters;2nativecancerrawprojections;full45traitQC/396replay incomplete'),
        (2,'PARTIAL','100lockedHEADversions available;100h2/1200rg printedlogs and603BH verified previously;nativefullchain incomplete'),
        (3,'PARTIAL','217classification retained23/18/17/159;0bothtraitindependent successes;MVPpreprocessing826027only'),
        (4,'METHOD_BLOCKED','No justified calibrated shared-sleep covariance;nominalheterogeneity remainsdiagnostic'),
        (5,'NO_GO_CURRENT_CONFIRMATORY_SLEEP','No eligible newprimaryquestion'),(6,'NO_EXPERIMENT_ADMITTED','0newcorrelations'),
        (7,'NOT_ADMITTED','Nolocalarchitectureexperiment'),(8,'NOT_ADMITTED','Nomolecularsignalormechansim'),
        (9,'PARTIAL','Frozen tables/sevenv1figures+onev3QCfigure;fullnativecomparisonsincomplete'),
        (10,'QUALIFIED_REVIEW_COMPLETE','Sevenv1separateroles+threev3roles;humanreviewunresolved'),
        (11,'NOT_READY','ScientificcompletenessandSLEEPnoveltyrequirementsunmet'),(12,'VERSIONED_EVIDENCE','V1/V2sealsunchanged;V3technicalevidence/tests/export')]])
    human='''# Human scientific review requirements v3

- Provide at least 3 GiB free internally, or a suitable compute host, before the unchanged full native runner is used. Source-only workload guards are separately justified and do not change this gate.
- Resolve the prostate historical container SHA649921... versus current public SHA89fa72... discrepancy. Original raw/provenance agree with the historical contract; source category counts materially differ. Current bytes cannot substitute for the locked source. No cause is established by matching HEAD or range checks.
- Confirm MVP confidence-interval level/construction, upstream effect allele semantics and EBI conversion; interpret P/CI discrepancies and R2>1 without post-hoc favorable filters. All supplied SE values are absent.
- Review clinical-insomnia phenotype equivalence, historical CIPHER implementation, MAC discrepancy and participant overlap. A preprocessed clinical transport source does not establish two-trait independent replication.
- Review high lipid cross-trait intercepts, MS/melanoma liability scale inputs and source prevalence/N assumptions, aligned jackknife covariance calibration, selection/winner's curse and the v1 novelty crosswalk. No heterogeneity or biological mechanism is confirmed.
- Choose any future primary estimand only after data, power, uncertainty and independence requirements are met. Current verdict remains NO-GO for a new confirmatory SLEEP expansion.

No manuscript, investigator contact, restricted access request, journal submission or main merge has been performed.
'''
    textfile(PACKAGE/'HUMAN_SCIENTIFIC_REVIEW_V3.md',human)
    final='''# Technical evidence handoff v3 — research only

Scientific status: incomplete and not ready for human manuscript authoring. This addendum preserves the separate 396/1,200/217 families and all sealed v1/v2 outputs. No new correlation or independent success is reported.

| Required item | Evidence disposition |
|---|---|
| A. Exact sources | Cumulative 45/45 raw and 42/43 container hashes available; v3 adds exact breast (4,764,280,363 bytes) and ovarian (3,774,309,689 bytes). Failed public prostate (2,534,497,125 bytes) preserved and rejected. All 100 immutable extension versions pass HEAD identity; bodies not newly verified. |
| B. Original artifacts | Prior 18/18 checkpoint files remain exact; cumulative 18/19 missing source-ledger entries resolved, 15 unique missing hashes recovered. Environment append-order discrepancy remains explicitly versioned. |
| C. Native reproduction | Prior insomnia–BMI control and three lipid rsID prefilters; v3 unchanged breast/ovarian raw projections reproduce every decompressed byte, line and CRC. MVP 19,703,815 rows → 826,027 HM3; source/build/allele/N and 64 samples independently checked. Full 45-trait QC/396, 100-trait QC/1,200 and 41-validation replay incomplete. |
| D. Differences | Cancer gzip serialization differs; decompressed raw contents identical. Independent review distinguishes header and deflate differences without inferring their cause. MVP 64/64 sampled Z match frozen binary-input arithmetic at .12g; exact decimal text differs at one final-digit boundary (maximum relative difference 4.88e-12). Prostate source SHA and category counts materially differ; no substitution. |
| E. Independent validation | 0 both-trait-independent successes. Historical 217 classification remains 23 qualified outcome-side positives, 18 concordant, 17 QC-ineligible, 159 without eligible external data; 41 estimates reused sleep GWAS. |
| F. New executed work | Exact-source transfers, two native raw projections, full MVP preprocessing, independent source/output streaming audit, 100-version HEAD checks, rejected-prostate metadata diagnostic and traceable QC figure. No new h2/rg, covariance-calibrated comparison, local or molecular experiment. |
| G. New scientific finding | No distinct new sleep–disease finding established; no new primary question admitted. |
| H. Prior findings | v1 supplement-based novelty crosswalk remains 15 substantially similar, 6 related, 2 unresolved among 23 positives. Morrison 2024, Goodman 2025, Sleep Chart 2026 and disease-specific prior evidence prevent first-ever claims. |
| I. Unresolved | Full native resource gate, historical prostate source, CI/effect semantics, h2/power, participant overlap, covariance calibration, intercept/liability issues and novelty. |
| J. Allocation | Keep one integrated research resource with separate families. Current evidence does not justify the intended SLEEP Original Article. After scientific validation, conditional specialist resource allocation such as SLEEP Advances is more defensible; a data-descriptor outlet requires suitable open/shareable data. |
| K. Human review | See HUMAN_SCIENTIFIC_REVIEW_V3.md plus all v1/v2 review requests. |
| L. Authoring clearance | NO. No manuscript text has been generated. |
| M. Artifacts | v3 scripts 19–27, tables/, figures/MVP_source_preprocessing_QC_v3.{png,pdf,svg}, logs/, manifests/, reviews/ and portable tests. Sealed v1/v2 remain siblings. Exact v3 Git commit/export proof is in the separate output audit receipt. |

Source status carries dated v1/v2 cryptographic proofs forward; new containers and projected raw inputs have current before/after receipts. Test-suite checks verify artifact contracts and never substitute for native scientific reproduction. The seven historical adversarial roles and three current methods/provenance/numerical roles are separately identified.

Full native preflight remains below 3 GiB internal free in three consecutive goal turns. Completed light work does not relax this requirement. The remaining native program requires resources and scientific review; goal completion is not claimed.
'''
    textfile(PACKAGE/'FINAL_EVIDENCE_HANDOFF_V3.md',final)
    textfile(PACKAGE/'README.md','''# Unified sleep research evidence v3

Research-only, separately versioned addendum to sealed v1/v2. Start with FINAL_EVIDENCE_HANDOFF_V3.md. Source/phenotype metadata and protected GWAS bodies remain on SSD; this package contains code, technical receipts, aggregate tables and one QC figure.

Two original source containers and raw projections pass; public prostate differs and is excluded. MVP materialization passes preprocessing-only review; all downstream inferential gates remain unchanged. Scientific goal incomplete; no manuscript-authoring clearance.

Run portable unittest discovery in tests/. Scientific executors require exact mounted SSD paths and frozen plans. Never rerun them into existing outputs. Use a new version and reviewed resource plan for any further execution. Full LDSC requires3GiB internal free. Reproduction tests/clean-checkout and Git export receipts are separate from source-to-result scientific reproduction.
''')
    save(PACKAGE/'logs/consolidation_receipt_v3.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),'evidence_state':state,
         'new_tables_sha256':{p.name:sha(p) for p in (PACKAGE/'tables').glob('*') if p.is_file()},'earlier_sealed_files_unchanged':278})
    print(json.dumps(state),flush=True)

if __name__=='__main__':main()
