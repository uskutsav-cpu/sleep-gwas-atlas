#!/usr/bin/env python3
"""Metadata/source-only binding comparison for canonical successor; no bodies."""
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4');R=P/'reviews'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def defs(p):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(p).read_text()).body if isinstance(n,ast.FunctionDef)}
def main():
 old=P/'manifests/native_canonical_calibration_plan_v4_3.json';new=P/'manifests/native_canonical_calibration_plan_v4_6.json'
 assert sha(old)=='cfe5ef24f3455cc0508785e69b291d023f524b2e4872754598b864365f884c21' and sha(new)=='3bbbaa75b5520e569c3386eb32cc28967cd5595ce1beb7dbea0cb7202c23b6ce'
 a=json.loads(old.read_text());b=json.loads(new.read_text());jobs=copy.deepcopy(a['jobs'])
 for j in jobs:
  j['output_dir']=j['output_dir'].replace(a['ssd_output_root']+'/',b['ssd_output_root']+'/',1)
  j['ldsc_args']=[x.replace(a['ssd_output_root']+'/',b['ssd_output_root']+'/',1) if x.startswith(a['ssd_output_root']+'/') else x for x in j['ldsc_args']]
 assert b['jobs']==jobs and len(jobs)==16
 changed=[k for k in set(a)|set(b) if a.get(k)!=b.get(k)]
 allowed={'jobs','self_path','ssd_output_root','dependencies_sha256','operational_revision','prepared_utc','preserved_v4_3_plan','preserved_v4_3_plan_sha256','global_resource_ledger_path','global_resource_ledger_sha256','global_reservation_bytes','terminal_protocol','resource_preflight'}
 assert set(changed)<=allowed,changed
 assert all(b['dependencies_sha256'].get(k)==v for k,v in a['dependencies_sha256'].items())
 text_equal={}
 for kind in ['capture','arithmetic']:
  x=P/'scripts'/('canonical_calibration_'+kind+'_v4_3.py');y=P/'scripts'/('canonical_calibration_'+kind+'_v4_6.py')
  text_equal[kind]=y.read_text().replace('_v4_6','_v4_3')==x.read_text();assert text_equal[kind]
 x=defs(P/'scripts/42_prepare_native_canonical_calibration_v4_3.py');y=defs(P/'scripts/42_prepare_native_canonical_calibration_v4_6.py')
 inherited=['monitor_helpers','snapshot','output_bytes','owned_group_exists','await_owned_cleanup','monitored_worker']
 assert all(x[k]==y[k] for k in inherited)
 commonold=(P/'scripts/canonical_calibration_common_v4_3.py').read_text();commonnew=(P/'scripts/canonical_calibration_common_v4_6.py').read_text()
 adjusted=commonnew.replace('_v4_6','_v4_3').replace('import native_stage_completion_v4_3 as stage_completion\n','')
 before="        p = PACKAGE / 'logs' / (stage + '_native_monitor_receipt_v4.json')\n        r = json.loads(p.read_text())"
 after='        p = stage_completion.stage_monitor_path(stage)\n        r = stage_completion.stage_monitor(stage)'
 assert adjusted.replace(after,before)==commonold
 draft=P/'manifests/native_canonical_calibration_plan_v4_4.json';assert sha(draft)=='9e21114163a68d39b9462ff80eb9ab01a0ca4bfd04c52acf670544d76cf46a32'
 helpersource=(P/'scripts/native_stage_completion_v4_3.py').read_text();assert "suffix = '_v4_3.json' if stage == 'extension' else '_v4.json'" in helpersource
 checked={};excluded={}
 for path,h in b['dependencies_sha256'].items():
  q=Path(path)
  if any(x in path for x in ['/ref/','/data/','coordinates_','coordinate_map','/native/']) or q.suffix not in ['.py','.json','.md','.tsv','.sha256','.sh']:
   excluded[path]='body/binary/native output outside narrow metadata review';continue
  assert q.stat().st_size<=12<<20,path
  assert sha(q)==h,path;checked[path]=h
 completion=P/'logs/native190_root_independent_completion_addendum_v4.json'
 assert sha(completion)=='afaf006dfc2a4ca7870558442c2ea478e102f6597e7eefbb47e95dc359d5722b'
 root_evidence=json.loads(completion.read_text())
 ledger=json.loads(Path(b['global_resource_ledger_path']).read_text());assert sha(b['global_resource_ledger_path'])==b['global_resource_ledger_sha256']=='04e2ea61d9a06281468674c1a8326367ffb7e507857814365e63108ee38ebbc7'
 assert sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes'] and ledger['reserved_total_bytes']+ledger['unallocated_margin_bytes']==300<<30
 seal=R/'independent_canonical_executor_prelaunch_review_v4_5_seal.json';assert sha(seal)=='9b4812b31a5f38b2eedd8876ae0162120c1d71427831b3ebf4d08d4880490982'
 assert all(sha(q)==h for q,h in json.loads(seal.read_text())['file_sha256'].items())
 current=(P/'scripts/42_prepare_native_canonical_calibration_v4_6.py').read_text();prior=(P/'scripts/42_prepare_native_canonical_calibration_v4_5.py').read_text()
 oldast=defs(P/'scripts/42_prepare_native_canonical_calibration_v4_5.py')['final_violations']
 import ast as astmod
 t=astmod.parse(current);v=next(n for n in t.body if isinstance(n,astmod.FunctionDef) and n.name=='final_violations');v.body=[n for n in v.body if not (isinstance(n,astmod.Assign) and any(isinstance(x,astmod.Name) and x.id=='campaign_bytes' for x in n.targets)) and not (isinstance(n,astmod.If) and isinstance(n.test,astmod.Compare) and isinstance(n.test.left,astmod.Name) and n.test.left.id=='campaign_bytes')]
 assert astmod.dump(v,include_attributes=False)==oldast
 assert b['arithmetic']==a['arithmetic'] and b['allow_41_covariance_outcomes'] is False and b['allow_calibrated_biological_p_values'] is False
 out=R/'independent_canonical_v4_6_binding_receipt.json'
 save(out,{'schema':'independent_canonical_v4_6_metadata_binding_check','checker_sha256':sha(__file__),'new_plan_sha256':sha(new),'old_plan_sha256':sha(old),'unexecuted_v4_4_plan_sha256':sha(draft),'all_binding_checks_pass':True,'native_control_jobs':16,'allowed_control_ids':b['allowed_control_ids'],'three_unique_source_hashes_unchanged':b['inputs_current_sha256']==a['inputs_current_sha256'],'exact_jobs_argv_except_private_paths':True,'capture_arithmetic_source_import_only_changes':text_equal,'inherited_monitor_functions_AST':inherited,'common_helper_changes_only_successor_paths_and_completion_selector':True,'top_level_plan_changes':sorted(changed),'arithmetic_tolerances_unchanged':b['arithmetic'],'verified_metadata_code_sha256':checked,'excluded_body_or_binary_sha256':excluded,'root190_completion_addendum_sha256':sha(completion),'root190_addendum_keys':sorted(root_evidence),'independent_real190_gate_executed':False,'GWAS_body_reads':0,'reference_body_reads':0,'actual_workers':0,'actual_fits':0,'actual_mutex_acquisitions':0,'execution_admission_granted':False,'preserved_v4_5_review_seal_sha256':sha(seal),'global_ledger_sha256':b['global_resource_ledger_sha256'],'global_reserved_bytes':ledger['reserved_total_bytes'],'global_margin_bytes':ledger['unallocated_margin_bytes'],'final_violations_only_adds_global300_meter':True})
 print(json.dumps({'receipt_sha256':sha(out),'metadata_files_verified':len(checked),'bodies_binaries_excluded':len(excluded)},indent=2))
if __name__=='__main__':main()
