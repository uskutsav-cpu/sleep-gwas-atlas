#!/usr/bin/env python3
"""Additive metadata map/mask/N-convention/CI proof after whole62 arithmetic."""
import csv,hashlib,json,math,os
from pathlib import Path
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4');R=P/'reviews';T=P/'tables';S=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
def sha(q):return hashlib.sha256(Path(q).read_bytes()).hexdigest()
def read(q):return json.loads(Path(q).read_text())
def tab(q):return list(csv.DictReader(Path(q).open(),delimiter='\t'))
def save(q,v):
 with Path(q).open('x') as f:json.dump(v,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def main():
 receipt=R/'independent_sensitivity_whole62_receipt_v4_6.json';assert sha(receipt)=='b758fa6f74db062873e597ce33d858a28080a67b42d0ff8b756f14b79c9f7ee2'
 proof=read(receipt);plan=read(S/'sensitivity_operational_plan_v4_6.json');original=P.parent/'sleep_unified_research_v1/tables/original_core_396.tsv';orig={(r['sleep_trait'],r['disease_trait']):r for r in tab(original)}
 oldrows=tab(T/'core_native_full_precision_rg_v4.tsv');bmap={(r['sleep_trait'],r['outcome_trait']):r for r in oldrows};rows=tab(T/'independent_sensitivity_62_full_precision_v4_6.tsv');checks=[]
 for a in plan['audit_jobs']:
  r=read(a['out_prefix']+'.intersection.json');assert r['dependency_sha256']==plan['dependencies_sha256'];checks.append(a['audit_id']+'/exact_full_dependency_map')
 for row in rows:
  j=next(j for j in plan['jobs'] if j['job_id']==row['job_id']);c=read(j['out_prefix']+'.full_precision.json');e=next(e for e in c['estimates'] if Path(e.get('p2',e.get('input'))).name.removesuffix('.sumstats.gz')==row['outcome_trait']);identity=next(x for x in c['final_intersections'] if Path(x.get('p2',x.get('input'))).name.removesuffix('.sumstats.gz')==row['outcome_trait'])
  objects=[('hsq1',e['hsq1']),('hsq2',e['hsq2']),('gencov',e['gencov'])] if j['kind']=='rg' else [('hsq',e)]
  for k,o in objects:
   assert o['n_blocks']==200 and o['constrain_intercept'] is False
   mask=identity.get('two_step_masks',{}).get(k) if identity.get('two_step_masks') is not None else identity.get('two_step_hsq_mask') if j['kind']=='h2' else None
   if mask is not None:assert o['twostep_filtered']==mask['total']-mask['count']
   else:assert o['twostep_filtered'] is None
  if j['kind']=='rg':
   base=bmap[(row['sleep_trait'],row['outcome_trait'])];o=orig[(row['sleep_trait'],row['outcome_trait'])]
   assert float(row['historical_family_q_value'])==float(o['fdr'])
   for k,x in [('historical_analysis_tier','analysis_tier'),('historical_interpretation_status','interpretation_status'),('historical_sleep_h2_verdict','sleep_h2_verdict'),('historical_outcome_h2_verdict','disease_h2_verdict'),('historical_exclusion_reason','disease_h2_qc_reason')]:assert row[k]==o[x]
   assert float(row['baseline_rg_ratio'])==float(base['rg'])
   assert float(row['baseline_rg_se'])==float(base['se'])
   for arm in ['baseline','sensitivity']:
    for k in ['hsq1','hsq2','gencov']:
     for f in ['tot','intercept']:
      point=float(row[arm+'_'+k+'_'+f]);se=float(row[arm+'_'+k+'_'+f+'_se']);row[arm+'_'+k+'_'+f+'_CI95_lower']=point-1.96*se;row[arm+'_'+k+'_'+f+'_CI95_upper']=point+1.96*se
   checks.append(row['sleep_trait']+'__'+row['outcome_trait']+'/frozen396_q_exclusions_exact')
 npath=P/'statistical_validation/processed_N_coding_diagnostics_v4.tsv';oldreceipt=P/'statistical_validation/qc_arithmetic_receipt_v4.json';old=read(oldreceipt);entry=next(x for x in old['outputs'] if x['path']==str(npath));assert sha(npath)==entry['sha256']=='412be00c7ce13fa0a0c31181cf430b046dd6b2f62122360305a25ec46e34b18f'
 n={r['trait']:r for r in tab(npath)};present=tab(T/'independent_sensitivity_liability_presentations_v4_6.tsv');diagnostics=[]
 for trait in ['ms','melanoma']:
  r=n[trait];assert r['source_sha256']==plan['input_sha256'][r['path']] and r['N_min']==r['N_max'] and r['distinct_N_count_capped_at_101']=='1' and r['gzip_read_to_EOF']=='True'
  sourceN=float(r['N_min']);base=next(x for x in present if x['trait']==trait and x['N_arm']=='original_effective_N' and x['P_convention']=='balanced_P_0_5');new=next(x for x in present if x['trait']==trait and x['N_arm']=='total_N_derivative' and x['P_convention']=='actual_original_case_fraction');cases=int(base['historical_case_count']);controls=int(base['historical_control_count']);total=cases+controls;Neff=4*cases*controls/total
  diagnostics.append({'trait':trait,'source_N_inventory_scope':'Previously sealed full processed-stream scan; not reread in this review','source_sha256':r['source_sha256'],'finite_source_N_constant':sourceN,'exact_count_based_Neff':Neff,'source_N_minus_exact_Neff':sourceN-Neff,'total_N':total,'N_scale_ratio_source_to_total':sourceN/total,'original_balanced_P_presented_h2':float(base['presented_h2']),'total_N_actual_P_presented_h2':float(new['presented_h2']),'observed_presented_h2_difference':float(new['presented_h2'])-float(base['presented_h2']),'observed_presented_relative_difference':float(new['presented_h2'])/float(base['presented_h2'])-1,'rounding_ratio_source_N_to_exact_Neff_minus1':sourceN/Neff-1,'no_tolerance_applied_to_cross_fit_difference':True,'variance_fraction_certified':False})
 output=T/'independent_sensitivity_62_full_precision_with_CI_v4_6.tsv';keys=list(dict.fromkeys(k for r in rows for k in r))
 with output.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);f.flush();os.fsync(f.fileno())
 out=R/'independent_sensitivity_whole62_addendum_receipt_v4_6.json';save(out,{'schema':'independent_whole62_identity_liability_addendum','checker_sha256':sha(__file__),'whole62_receipt_sha256':sha(receipt),'full52_audit_dependency_maps_exact':True,'native200_block_free_intercepts_and_mask_counts_match':True,'historical60_frozen396_q_and_exclusions_exact':True,'baseline60_fullprecision_table_values_match':True,'check_count':len(checks),'checks':checks,'original396_table_sha256':sha(original),'processed_N_inventory_sha256':sha(npath),'original_qc_arithmetic_receipt_sha256':sha(oldreceipt),'N_and_liability_descriptive_diagnostics':diagnostics,'enriched62_table_sha256':sha(output),'source_GWAS_reference_bodies_reread':False,'actual_mutex_acquisitions':0,'actual_workers':0,'actual_fits':0,'cross_fit_difference_P_values':0,'new_QC_admission':False})
 print(json.dumps({'receipt_sha256':sha(out),'enriched_table_sha256':sha(output),'diagnostics':diagnostics},indent=2))
if __name__=='__main__':main()
