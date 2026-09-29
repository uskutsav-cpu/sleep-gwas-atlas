#!/usr/bin/env python3
"""Validate the Brain6-only, outcome-blinded sleep GWAS power screen."""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'brain6/results/power_optimized_sensitivity_v1'
def sha(path: Path)->str: return hashlib.sha256(path.read_bytes()).hexdigest()
def read_tsv(path: Path):
 with path.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))
def validate(root: Path=ROOT)->dict:
 base=root/'brain6/results/power_optimized_sensitivity_v1'; prov=json.loads((base/'provenance.json').read_text())
 if prov['scope']!='Brain6 only; outcome-blinded candidate screen; no frailty/FI results or sensitivity receipts included': raise ValueError('scope mismatch')
 if prov['builder_sha256']!=sha(root/prov['builder_script']): raise ValueError('builder hash mismatch')
 for p,h in prov['inputs'].items():
  if sha(root/p)!=h: raise ValueError(f'input hash mismatch: {p}')
 for p,h in prov['outputs'].items():
  if sha(base/p)!=h: raise ValueError(f'output hash mismatch: {p}')
 traits=read_tsv(base/'current_sleep_trait_audit.tsv'); cand=read_tsv(base/'candidate_replacement_screen.tsv')
 if {r['trait'] for r in traits}!={'insomnia','longsleep'}: raise ValueError('must audit only the two sleep GWAS in the seven-input family')
 expected={'insomnia':(1824,671,651),'longsleep':(1204,1291,1271)}
 for r in traits:
  if tuple(int(r[k]) for k in ('processable_tested','not_run','low_local_h2_not_run'))!=expected[r['trait']]: raise ValueError(f"canonical receipt counts changed for {r['trait']}")
 if int(prov['decision_rule']['minimum_additional_loci'])!=250 or prov['decision_rule']['material_improvement_strict_gate_eligible_loci_percentage_points']!=10: raise ValueError('promotion rule changed')
 if len(cand)<3 or prov['status']!='COMPLETE_SCREEN_NO_PROMOTION': raise ValueError('candidate screen incomplete or status mismatch')
 if prov.get('trait_only_lava_screen')!='COMPLETE_2495_LOCI_NO_PROMOTION': raise ValueError('candidate trait-only screen is not recorded as complete')
 screen=base/'lava_trait_screen_v1'
 aggregate_summary=json.loads((screen/'summary.json').read_text())
 aggregate_provenance=json.loads((screen/'aggregate.provenance.json').read_text())
 completion_provenance=json.loads((base/'completed_trait_only_screen_comparison_v1.provenance.json').read_text())
 aggregate_rows=read_tsv(screen/'sleep_duration_continuous_lava_univariate.tsv')
 if (aggregate_summary.get('status')!='COMPLETE_SCREEN_NO_FAILED_ROWS' or aggregate_summary.get('rows')!=2495 or
     aggregate_summary.get('tested')!=1619 or aggregate_summary.get('not_run')!=876 or aggregate_summary.get('failed')!=0 or
     aggregate_summary.get('strict_gate_pass')!=7 or aggregate_summary.get('canonical_longsleep_strict_gate_pass')!=1 or
     aggregate_summary.get('additional_strict_gate_pass')!=6 or aggregate_summary.get('minimum_additional_for_promotion')!=250 or
     aggregate_summary.get('eligible_for_full_sensitivity_screen_gate') is not False or len(aggregate_rows)!=2495):
  raise ValueError('candidate trait-only screen summary is incomplete or would permit unqualified promotion')
 if (aggregate_provenance.get('status')!='COMPLETE_SCREEN_NO_FAILED_ROWS' or
     aggregate_provenance.get('base_config_sha256')!=sha(screen/'config.json') or
     aggregate_provenance.get('outputs',{}).get('brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/sleep_duration_continuous_lava_univariate.tsv')!=sha(screen/'sleep_duration_continuous_lava_univariate.tsv') or
     completion_provenance.get('status')!='COMPLETE_SCREEN_NO_PROMOTION' or
     completion_provenance.get('summary',{}).get('full_sensitivity_eligible') is not False or
     completion_provenance.get('decision_rule',{}).get('downstream_association_results_consulted') is not False):
  raise ValueError('candidate screen provenance or outcome-blind promotion decision is invalid')
 complete=read_tsv(base/'current_sleep_trait_audit_complete.tsv')
 overlap_prov=json.loads((base/'canonical_sleep_reference_overlap_v1.provenance.json').read_text())
 overlap=read_tsv(base/'canonical_sleep_reference_overlap_v1.tsv')
 if (overlap_prov.get('status')!='PASS_HASH_VERIFIED_EXACT_ID_OVERLAP' or
     overlap_prov.get('builder_sha256')!=sha(root/overlap_prov.get('builder_path','')) or
     any(not Path(path).is_file() or sha(Path(path))!=digest for path,digest in overlap_prov.get('inputs_sha256',{}).items()) or
     any(sha(base/name)!=digest for name,digest in overlap_prov.get('outputs_sha256',{}).items()) or
     len(complete)!=2 or len(overlap)!=2):
  raise ValueError('canonical sleep reference-overlap audit is incomplete or hash-invalid')
 expected_overlap={'insomnia':(6077635,6061267),'longsleep':(6549769,6057783)}
 for row in complete:
  source_n,overlap_n=expected_overlap[row['trait']]
  if (int(row['usable_variant_count'])!=source_n or int(row['frozen_reference_overlap_count'])!=overlap_n or
      int(row['frozen_reference_nonoverlap_count'])!=source_n-overlap_n or
      abs(float(row['frozen_reference_overlap_pct'])-100*overlap_n/source_n)>1e-6):
   raise ValueError(f"source/reference overlap changed for {row['trait']}")
 if {row['trait'] for row in overlap}!={'insomnia','longsleep'}:
  raise ValueError('reference-overlap table does not cover both Brain6 sleep GWAS')
 corrected_path=base/'current_sleep_trait_audit_complete_v2.tsv'
 corrected=read_tsv(corrected_path)
 corrected_prov=json.loads((base/'current_sleep_trait_audit_complete_v2.provenance.json').read_text())
 canonical=read_tsv(base/'canonical_family_trait_power_audit_v1.tsv')
 canonical_by_trait={row['trait_id']:row for row in canonical}
 if (corrected_prov.get('status')!='PASS_CORRECTED_CANONICAL_METRICS' or
     corrected_prov.get('builder_sha256')!=sha(root/corrected_prov.get('builder_path','')) or
     any(not (root/path).is_file() or sha(root/path)!=digest for path,digest in corrected_prov.get('inputs_sha256',{}).items()) or
     any(sha(base/name)!=digest for name,digest in corrected_prov.get('outputs_sha256',{}).items()) or
     len(corrected)!=2):
  raise ValueError('corrected canonical sleep metrics are incomplete or hash-invalid')
 for row in corrected:
  source=canonical_by_trait[row['trait']]
  if (row['snp_h2'].split()[0]!=source['SNP_h2_reported_in_locked_global_map'] or
      row['local_h2_median_tested']!=source['local_h2_median_tested'] or
      row['local_h2_q25_tested']!=source['local_h2_q25_tested'] or
      row['local_h2_q75_tested']!=source['local_h2_q75_tested'] or
      row['local_h2_quartile_method']!='Linear interpolation; TESTED loci only'):
   raise ValueError(f"corrected metrics differ from canonical audit for {row['trait']}")
 candidates_v2_path=base/'candidate_replacement_screen_complete_v2.tsv'
 candidates_v2=read_tsv(candidates_v2_path)
 candidates_v2_prov=json.loads((base/'candidate_replacement_screen_complete_v2.provenance.json').read_text())
 v2_by_name={row['candidate']:row for row in candidates_v2}
 if (candidates_v2_prov.get('status')!='PASS_COMPLETE_SIX_CANDIDATE_SCREEN' or
     candidates_v2_prov.get('builder_sha256')!=sha(root/candidates_v2_prov.get('builder_path','')) or
     any(not (root/path).is_file() or sha(root/path)!=digest for path,digest in candidates_v2_prov.get('inputs_sha256',{}).items()) or
     any(sha(base/name)!=digest for name,digest in candidates_v2_prov.get('outputs_sha256',{}).items()) or
     len(candidates_v2)!=6 or
     candidates_v2_prov.get('full_sensitivity_candidates_eligible')!=0):
  raise ValueError('complete candidate screen is incomplete or hash-invalid')
 if not v2_by_name['Dashti 2019 continuous self-reported sleep duration']['decision'].startswith('NO FULL-SENSITIVITY PROMOTION'):
  raise ValueError('completed continuous-duration screen decision is stale')
 if not v2_by_name['Portas et al. 2026 device-measured sleep GWAS']['decision'].startswith('EXCLUDE from the long-sleep power screen'):
  raise ValueError('device-derived sleep candidate was not excluded on phenotype/power grounds')
 decision=json.loads((root/'work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json').read_text())
 if decision['overall_status']!='FAILED_QC_NOT_PROMOTED' or decision['verified_cells']!=17465 or decision['planned_loci']!=2495: raise ValueError('canonical invariant changed')
 mat=json.loads((base/'lava_inputs_v1/materialization.provenance.json').read_text()); rec=json.loads((base/'normalized/sleep_duration_continuous_dashti_2019/receipt.json').read_text())
 if rec['status']!='COMPLETE' or rec['scientific_status']!='PASS' or mat['loci']!=2495 or mat['frozen_reference_unique_overlap']!=6093758: raise ValueError('candidate input/coverage gate changed')
 return {'status':'PASS_BRAIN6_POWER_AUDIT','sleep_traits':len(traits),'candidates_screened':len(candidates_v2),'candidate_loci_materialized':mat['loci'],'candidate_reference_overlap_variants':mat['frozen_reference_unique_overlap'],'canonical_sleep_reference_overlap':{row['trait']:int(row['frozen_reference_overlap_count']) for row in corrected},'candidate_screen':{'tested':aggregate_summary['tested'],'not_run':aggregate_summary['not_run'],'strict_gate_pass':aggregate_summary['strict_gate_pass'],'full_sensitivity_promoted':False},'canonical_cells':decision['verified_cells'],'canonical_decision':decision['overall_status'],'corrected_sleep_metric_rows':len(corrected)}
if __name__=='__main__': print(json.dumps(validate(),sort_keys=True))
