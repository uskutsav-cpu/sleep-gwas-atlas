#!/usr/bin/env python3
"""Independent stdlib captured/delete-centered review; no scientific worker/import."""
import csv,hashlib,json,math,os,resource,statistics,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4');R=P/'reviews';T=P/'tables';ROOT=P.parent
S=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
PLAN=S/'sensitivity_operational_plan_v4_6.json';PIN='b8f9c4c0dfecf4f5f4559df3f468b1f2128ef932caece54691124a8788432efb'
STAGE=S/'receipts_v4/sensitivity_execution_receipt_v4_6.json';STAGEPIN='c674c83b7d6f6dc25bbaa4b6b5388bfeade67a23f20f9bb3df50f9a67d920460'
SEAL=S/'receipts_v4/stage_terminal_seal_v4_6.json';SEALPIN='054e3b41fa26729be23d6420544f14e4b5e303ae8bd20b2287e81c0f76352669'
START=time.monotonic();CHECKS=[];HASHES={};MAXERR={};WARN=[];BODY_READS=0
BUDGET=json.loads((R/'independent_sensitivity_whole62_resource_plan_v4_6.json').read_text())
def usage():return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
def bounded(q):
 q=Path(q);assert q.stat().st_size<=BUDGET['single_read_file_limit_bytes'],q
 assert usage()<BUDGET['RSS_ceiling_bytes'] and time.monotonic()-START<BUDGET['per_process_seconds_ceiling']
 return q

def sha(q):
 q=bounded(q);h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def regular(q,wanted=None):
 q=Path(q);assert q.is_file() and not q.is_symlink() and not any(p.is_symlink() for p in q.parents),q
 actual=sha(q)
 if wanted:assert actual==wanted,(q,actual,wanted)
 HASHES[str(q)]=actual;return actual
def read(q):return json.loads(bounded(q).read_text())
def save(q,v):
 with Path(q).open('x') as f:json.dump(v,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def tsv(q,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(q).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v,sort_keys=True,separators=(',',':')) if isinstance(v,(dict,list)) else ('NULL' if v is None else v) for k,v in r.items()})
  f.flush();os.fsync(f.fileno())
def check(name,good,detail=None):
 CHECKS.append({'name':name,'pass':bool(good),'detail':detail});assert good,(name,detail)
def close(name,a,b,p=False,presentation=False):
 if not all(isinstance(x,(float,int)) and math.isfinite(x) for x in [a,b]):check(name,False,[a,b])
 diff=abs(a-b);category=name.rsplit('/',1)[-1];MAXERR[category]=max(MAXERR.get(category,0),diff)
 check(name,math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12 if presentation else (1e-300 if p else 1e-15)),{'absolute_error':diff});return diff

def vector(v):
 x=[row[0] if isinstance(row,list) and len(row)==1 else row for row in v]
 assert len(x)==200 and all(isinstance(t,(int,float)) and math.isfinite(t) for t in x)
 return x

def variance(v):
 center=math.fsum(v)/200
 return 199/200*math.fsum((x-center)**2 for x in v)
def se(v):return math.sqrt(variance(v))
def object_checks(obj,label):
 d=vector(obj['tot_delete_values']);i=vector(obj['intercept_delete_values'])
 close(label+'/total_se',se(d),obj['tot_se']);close(label+'/total_cov',variance(d),obj['tot_cov']);close(label+'/intercept_se',se(i),obj['intercept_se'])
 jk=obj['jknife'];de=jk['delete_values'];assert len(de)==200 and all(len(x)==2 for x in de)
 assert i==[x[1] for x in de]
 centers=[math.fsum(x[k] for x in de)/200 for k in range(2)]
 for k in range(2):
  values=[x[k] for x in de];close(label+'/jackknife_se_'+str(k),se(values),jk['jknife_se'][0][k]);close(label+'/jackknife_mean_'+str(k),math.fsum(200*jk['est'][0][k]-199*x for x in values)/200,jk['jknife_est'][0][k])
  for l in range(2):close(label+'/jackknife_cov_'+str(k)+str(l),199/200*math.fsum((x[k]-centers[k])*(x[l]-centers[l]) for x in de),jk['jknife_cov'][k][l])
 return d

def identity_fields(kind):return ['final_ordered_SNP_count','final_ordered_SNP_sha256']+(['final_ordered_SNP_allele_compatibility_sha256','final_ordered_aligned_Z1_Z2_float64_big_endian_sha256'] if kind=='rg' else ['final_ordered_Z_float64_big_endian_sha256'])
def outcome(path):return Path(path).name.removesuffix('.sumstats.gz')
def tab(q):regular(q);return list(csv.DictReader(bounded(q).open(),delimiter='\t'))
def arguments_match(args,argv):
 if isinstance(args,list):assert args==argv;return
 i=0
 while i<len(argv):
  option=argv[i][2:].replace('-','_');actual=args[option]
  if isinstance(actual,bool):assert actual is True,(option,actual);i+=1;continue
  value=argv[i+1]
  if isinstance(actual,(int,float)):assert float(value)==actual,(option,value,actual)
  else:assert str(actual)==value,(option,value,actual)
  i+=2
 assert args['no_intercept'] is False and args['no_check_alleles'] is False and args['print_delete_vals'] is True
 assert args['intercept_h2'] is None and args['intercept_gencov'] is None and args['n_blocks']==200

def main():
 regular(PLAN,PIN);plan=read(PLAN);regular(STAGE,STAGEPIN);stage=read(STAGE);regular(SEAL,SEALPIN);seal=read(SEAL)
 pending=S/'receipts_v4/stage_pending_v4_6.json'
 check('terminal_pending_and_addenda_absent',not pending.exists() and not pending.is_symlink() and not Path(str(SEAL)+'.failure.json').exists() and not Path(str(SEAL)+'.failure.json').is_symlink() and not Path(str(STAGE)+'.failure.json').exists() and not Path(str(STAGE)+'.failure.json').is_symlink())
 check('exact_terminal_binding',seal['status']=='REVIEWED_STAGE_TERMINAL_SEAL' and seal['binding']=={'plan_sha256':PIN,'executor_sha256':sha(P/'scripts/sensitivity_executor_v4_6.py'),'merge_only':False} and seal['result_receipt_sha256']=={str(STAGE):STAGEPIN} and seal['pending_path']==str(pending))
 check('stage_exact_family',stage['plan_sha256']==PIN and len(stage['completed_jobs'])==26 and stage['new_fit_count']==62 and stage['merge_only'] is False and stage['historical_protocol_and_threshold_changes']==0 and len(plan['jobs'])==26 and sum(j['estimates'] for j in plan['jobs'])==62 and len(plan['audit_jobs'])==52)
 check('602_artifact_inventory',len(stage['consumed_output_sha256'])==602 and sum(Path(p).stat().st_size for p in stage['consumed_output_sha256'])<BUDGET['expected_consumed_total_bytes_upper_bound'])
 for q,h in stage['consumed_output_sha256'].items():regular(q,h)
 clock=S/'receipts_v4/stage_clock_v4_6.json';clockh=regular(clock,stage['stage_clock_sha256']);check('clock_immutable_consumed',stage['consumed_output_sha256'][str(clock)]==clockh and read(clock)['plan_sha256']==PIN)
 for q,h in stage['baseline_gate_receipt_sha256'].items():regular(q,h)
 for q,h in plan['dependencies_sha256'].items():
  if any(x in q for x in ['/ref/','/data/','/native/']) or Path(q).suffix not in ['.py','.md','.json','.tsv','.sha256','.sh']:continue
  regular(q,h)
 protocol=P/'FROZEN_ESTIMATOR_SENSITIVITY_PROTOCOL_v1.md';regular(protocol,'ba7fb9464054c79a38625907544e9871cb6a5ee12c5bb46f5f116692a4061a79')
 manifest=read(plan['scientific_member_manifest']);regular(plan['scientific_member_manifest'],plan['scientific_member_manifest_sha256']);check('frozen62_matrix',[{k:j[k] for k in ['job_id','kind','inputs','options','estimates']} for j in plan['jobs']]==manifest['jobs'])
 proof=read(S/'proofs/final_stock_intersection_comparison_v4_6.json');regular(S/'proofs/final_stock_intersection_comparison_v4_6.json',stage['final_intersection_proof_sha256']);check('intersection_receipt',proof['plan_sha256']==PIN and proof['paired_estimates']==62 and proof['identity_instances']==124 and proof['estimator_calls']==0)
 deriv=read(plan['derivative_proof']);regular(plan['derivative_proof'],plan['derivative_proof_sha256']);check('full_stream_derivative_proof',deriv['status']=='RESULT_FREE_ONLY_FINITE_N_CHANGED' and [x['new_finite_N'] for x in deriv['results']]==[376169,290130] and all(x['original_and_derivative_invariant_match'] and x['literal_missingness_preserved'] and x['source_sha256_before']==x['source_sha256_after'] and x['full_stream_verified_rows']==1217311 and plan['input_sha256'][x['source']]==x['source_sha256_before'] and plan['input_sha256'][x['derivative']]==x['derivative_sha256'] for x in deriv['results']))
 baseline_rg=tab(T/'core_native_full_precision_rg_v4.tsv');baseline_h=tab(T/'core_native_full_precision_h2_v4.tsv');by_pair={(x['sleep_trait'],x['outcome_trait']):x for x in baseline_rg};by_trait={x['trait_id']:x for x in baseline_h}
 check('baseline_tables_cardinality',len(baseline_rg)==396 and len(baseline_h)==45)
 panel=tab(ROOT/'config/analysis_panel.tsv');casecounts={x['trait_id']:(int(x['ncase']),int(x['ncontrol'])) for x in panel if x['trait_id'] in ['ms','melanoma']}
 audits={};workers={};command_order=[]
 for a in plan['audit_jobs']:
  wp=S/'receipts_v4'/(a['audit_id']+'.worker.json');w=read(wp);workers[a['audit_id']]=w
  expected=[plan['python'],'-B',str(P/'scripts/37_stock_ldsc_intersection_audit.py'),'--plan',str(PLAN),'--plan-sha',PIN,'--audit-id',a['audit_id']]
  check(a['audit_id']+'/exact_worker_command',w['command']==expected)
  r=read(a['out_prefix']+'.intersection.json');check(a['audit_id']+'/stock_audit_receipt',r['audit_id']==a['audit_id'] and r['arm']==a['arm'] and r['plan_sha256']==PIN and r['estimator_calls']==0 and r['estimator_constructors_forbidden'] is True and len(r['final_intersections'])==a['expected_identities'] and r['source_sha256_before']==r['source_sha256_after']=={p:plan['input_sha256'][p] for p in a['inputs']})
  arguments_match(r['arguments'],a['ldsc_args']);audits[a['audit_id']]=r
 for j in plan['jobs']:
  w=read(S/'receipts_v4'/(j['job_id']+'.worker.json'));workers[j['job_id']]=w
  expected=[plan['python'],'-u','-B',str(P/'scripts/38_sensitivity_ldsc_capture.py'),'--plan',str(PLAN),'--plan-sha',PIN,'--job-id',j['job_id']]+j['ldsc_args']
  check(j['job_id']+'/exact_fit_command',w['command']==expected)
 for key,w in workers.items():
  check(key+'/worker_complete',w['status']=='WORKER_COMPLETE_VERIFIED' and w['plan_sha256']==PIN and w['returncode']==0 and w['stop_reason'] is None and w['plan_unchanged'] is True and w['owned_cleanup_verified'] is True and w['process_group_teardown']['remaining_group_members']==[] and w['metadata_errors']==[])
  for q,h in w['output_sha256'].items():check(key+'/consumed_worker_output',stage['consumed_output_sha256'].get(q)==h)
 check('all_audits_completed_before_first_fit',max(workers[a['audit_id']]['completed_utc'] for a in plan['audit_jobs'])<min(workers[j['job_id']]['resource_preflight']['recorded_utc'] for j in plan['jobs']))
 check('78_exact_worker_members',{str(S/'receipts_v4'/(k+'.worker.json')) for k in workers}=={q for q in stage['consumed_output_sha256'] if q.endswith('.worker.json')})
 rows=[];arithmetic=[];pairs=[];standalone={};stock_arrays=0;baseline_capture_hashes={}
 for j in plan['jobs']:
  cp=Path(j['out_prefix']+'.full_precision.json');c=read(cp);captured=next(x for x in stage['completed_jobs'] if x['job_id']==j['job_id'])
  check(j['job_id']+'/capture_identity',c['job_id']==j['job_id'] and c['plan_sha256']==PIN and len(c['estimates'])==j['estimates'] and captured['estimates']==j['estimates'] and captured['full_precision_sha256']==HASHES[str(cp)] and c['input_sha256_before']==c['input_sha256_after']=={q:plan['input_sha256'][q] for q in j['inputs']} and c['libraries']=={'numpy':'1.21.5','pandas':'1.3.3','scipy':'1.7.3'} and c['python'].startswith('3.9.23'))
  arguments_match(c['arguments'],j['ldsc_args']);log=Path(j['out_prefix']+'.log');regular(log,c['stock_log_sha256'])
  b_a=next(a for a in plan['audit_jobs'] if a['sensitivity_job_id']==j['job_id'] and a['arm']=='baseline');s_a=next(a for a in plan['audit_jobs'] if a['sensitivity_job_id']==j['job_id'] and a['arm']=='sensitivity')
  bb=audits[b_a['audit_id']]['final_intersections'];ss=audits[s_a['audit_id']]['final_intersections'];check(j['job_id']+'/native_final_intersection_matches_stock_audit',c['final_intersections']==ss)
  for path,h in c['stock_delete_array_sha256'].items():regular(path,h)
  for index,(e,x,y) in enumerate(zip(c['estimates'],bb,ss)):
   label=j['job_id']+'/'+str(index);fields=identity_fields(j['kind']);sameN=x['final_ordered_N1_N2_float64_big_endian_sha256' if j['kind']=='rg' else 'final_ordered_N_float64_big_endian_sha256']==y['final_ordered_N1_N2_float64_big_endian_sha256' if j['kind']=='rg' else 'final_ordered_N_float64_big_endian_sha256']
   check(label+'/exact_paired_identity',all(x[k]==y[k] for k in fields) and sameN==j['job_id'].startswith('lipid_two_step_'))
   assert x.get('p2',x.get('input')).split('/')[-1]==y.get('p2',y.get('input')).split('/')[-1]
   pair=next(p for p in proof['pairs'] if p['job_id']==j['job_id'] and p['outcome']==Path(y.get('p2',y.get('input'))).name);check(label+'/retained_pair_proof',pair['baseline']==x and pair['sensitivity']==y and pair['N_identity_same']==sameN)
   objects=[e['hsq1'],e['hsq2'],e['gencov']] if j['kind']=='rg' else [e]
   deletes=[object_checks(o,label+'/'+str(i)) for i,o in enumerate(objects)]
   # Independent text stock .delete vectors must reproduce captured totals exactly.
   for role,vals in zip(['hsq1','hsq2','gencov'] if j['kind']=='rg' else ['h2'],deletes):
    if j['kind']=='rg':suffix=Path(e['p1']).name+'_'+Path(e['p2']).name+'.'+role+'.delete';path=next(q for q in c['stock_delete_array_sha256'] if q.endswith(suffix))
    else:path=j['out_prefix']+'.delete'
    stock=[float(line) for line in Path(path).read_text().splitlines() if line.strip()];check(label+'/'+role+'/stock_tot_delete_exact',stock==vals);stock_arrays+=1
   if j['kind']=='h2':
    path=j['out_prefix']+'.part_delete';stock=[[float(v) for v in line.split()] for line in Path(path).read_text().splitlines() if line.strip()];check(label+'/stock_part_delete_exact',stock==e['part_delete_values']);stock_arrays+=1
   o=outcome(e.get('p2',e.get('input')));sleep=outcome(e['p1']) if j['kind']=='rg' else ''
   base=by_pair[(sleep,o)] if j['kind']=='rg' else by_trait[o];baseline_path=base['native_full_precision_path'];baseline_sha=base['native_full_precision_sha256'];regular(baseline_path,baseline_sha);baseline_capture_hashes[baseline_path]=baseline_sha
   bc=read(baseline_path);baseline_exec_path=baseline_path.replace('.full_precision.json','.execution_receipt.json');baseline_execution=read(baseline_exec_path);regular(baseline_exec_path,stage['baseline_gate_receipt_sha256'][baseline_exec_path]);check(label+'/baseline_source_receipt_maps',baseline_execution['input_sha256']==baseline_execution['input_sha256_after'] and baseline_execution['returncode']==0 and baseline_execution['all_output_sha256'][baseline_path]==baseline_sha)
   baseline_log=bc['arguments']['out']+'.log';regular(baseline_log,baseline_execution['all_output_sha256'][baseline_log]);baseline_warnings=[line for line in Path(baseline_log).read_text().splitlines() if any(t in line.upper() for t in ['WARNING','TRACEBACK','ERROR'])]
   be=next(q for q in bc['estimates'] if outcome(q.get('p2',q.get('input')))==o)
   r={'job_id':j['job_id'],'kind':j['kind'],'sleep_trait':sleep,'outcome_trait':o,'scope':'DESCRIPTIVE_PREDECLARED_SENSITIVITY_ONLY','status':e.get('status','NATIVE_H2_ESTIMATE_RETURNED'),'baseline_capture_path':baseline_path,'baseline_capture_sha256':baseline_sha,'sensitivity_capture_path':str(cp),'sensitivity_capture_sha256':HASHES[str(cp)],'baseline_input_sha256':{q:baseline_execution['input_sha256'][q] for q in [be.get('p1'),be.get('p2',be.get('input'))] if q is not None},'baseline_execution_receipt_path':baseline_exec_path,'baseline_execution_receipt_sha256':HASHES[baseline_exec_path],'sensitivity_input_sha256':c['input_sha256_before'],'source_hash_verification_scope':'Stored native before/after maps and previously sealed full-stream derivative proof; no GWAS/reference bodies reread by reviewer','baseline_options':bc['arguments'],'sensitivity_options':c['arguments'],'baseline_warnings':baseline_warnings,'baseline_warning_scope':'Full retained original native batch; pair-specific attribution not inferred','sensitivity_warnings':c['warnings_and_errors'],'final_SNP_count':y['final_ordered_SNP_count'],'final_SNP_sha256':y['final_ordered_SNP_sha256'],'allele_compatibility_sha256':y.get('final_ordered_SNP_allele_compatibility_sha256'),'baseline_audit_sha256':HASHES[b_a['out_prefix']+'.intersection.json'],'sensitivity_audit_sha256':HASHES[s_a['out_prefix']+'.intersection.json'],'baseline_audit_is_fit_instrumentation':False,'N_hash_equal':sameN,'baseline_effective_two_step':x['effective_two_step'],'sensitivity_effective_two_step':y['effective_two_step'],'sensitivity_two_step_masks':y.get('two_step_masks',y.get('two_step_hsq_mask')),'final_Z_filter_removed':y['final_Z_filter_removed'],'historical_family_q_value':base.get('frozen_fdr_original_complete_family'),'historical_primary_phase1_q_value':base.get('original_fdr_primary_phase1'),'historical_analysis_tier':base.get('original_analysis_tier'),'historical_interpretation_status':base.get('original_interpretation_status'),'historical_sleep_h2_verdict':base.get('original_sleep_h2_verdict'),'historical_outcome_h2_verdict':base.get('original_disease_h2_verdict',base.get('original_verdict')),'historical_exclusion_reason':base.get('original_disease_h2_qc_reason',base.get('original_qc_reason')),'new_BH':False,'difference_P_value':None,'cross_fit_delete_alignment_certified':False}
   if j['kind']=='rg':
    check(label+'/returned_finite',e['status']=='NATIVE_ESTIMATE_RETURNED' and all(math.isfinite(e[k]) for k in ['rg_ratio','rg_se','rg_jknife','z','p']))
    a,b,g=objects;ratio=g['tot']/math.sqrt(a['tot']*b['tot']);rd=[g/math.sqrt(a*b) for a,b,g in zip(*deletes)]
    for k,v in [('rg_ratio',ratio),('rg_se',se(rd)),('rg_jknife',math.fsum(200*ratio-199*x for x in rd)/200),('z',ratio/se(rd)),('p',math.erfc(abs(ratio/se(rd))/math.sqrt(2)))]:close(label+'/'+k,v,e[k],p=k=='p')
    for arm,item in [('baseline',be),('sensitivity',e)]:
     for k in ['rg_ratio','rg_jknife','rg_se','z','p']:r[arm+'_'+k]=item[k]
     r[arm+'_CI95_lower']=item['rg_ratio']-1.96*item['rg_se'];r[arm+'_CI95_upper']=item['rg_ratio']+1.96*item['rg_se']
     for name in ['hsq1','hsq2','gencov']:
      for k in ['tot','tot_se','intercept','intercept_se','mean_chisq','mean_z1z2','lambda_gc','ratio','ratio_se']:r[arm+'_'+name+'_'+k]=item[name].get(k)
      r[arm+'_'+name+'_h2_Z']=item[name]['tot']/item[name]['tot_se'] if name!='gencov' else None
     r[arm+'_pairwise_h2_Z_ge4_both']=all(item[n]['tot']/item[n]['tot_se']>=4 for n in ['hsq1','hsq2']);r[arm+'_pairwise_h2_intercept_le1_2_both']=all(item[n]['intercept']<=1.2 for n in ['hsq1','hsq2'])
    r['delta_rg']=e['rg_ratio']-be['rg_ratio'];r['delta_rg_SE']=e['rg_se']-be['rg_se'];r['intercept_boundary_changed']=r['baseline_pairwise_h2_intercept_le1_2_both']!=r['sensitivity_pairwise_h2_intercept_le1_2_both'];r['h2_Z_boundary_changed']=r['baseline_pairwise_h2_Z_ge4_both']!=r['sensitivity_pairwise_h2_Z_ge4_both']
   else:
    standalone[o]=e
    for arm,item in [('baseline',be),('sensitivity',e)]:
     for k in ['tot','tot_se','intercept','intercept_se','mean_chisq','lambda_gc','ratio','ratio_se']:r[arm+'_'+k]=item[k]
     r[arm+'_h2_Z']=item['tot']/item['tot_se'];r[arm+'_h2_Z_ge4']=r[arm+'_h2_Z']>=4;r[arm+'_intercept_le1_2']=item['intercept']<=1.2;r[arm+'_CI95_lower']=item['tot']-1.96*item['tot_se'];r[arm+'_CI95_upper']=item['tot']+1.96*item['tot_se']
   rows.append(r);arithmetic.append({'identity':sleep+'__'+o if sleep else o,'job_id':j['job_id'],'all_within_fit_checks_pass':True})
 check('all62_rows_retained',len(rows)==62 and sum(r['kind']=='rg' for r in rows)==60 and len(standalone)==2 and stock_arrays==184)
 presentations=[]
 for trait,new in standalone.items():
  old=by_trait[trait];K=float(old['population_prevalence_argument']);oldP=float(old['sample_prevalence_argument']);cases,controls=casecounts[trait];actualP=cases/(cases+controls);assert float(next(j for j in plan['jobs'] if j['job_id']=='binary_total_N_h2_'+trait)['options'][1])==actualP
  raw=float(old['h2_observed_full_precision']);rawse=float(old['h2_observed_se_full_precision']);nd=statistics.NormalDist();z=nd.inv_cdf(1-K);density=math.exp(-z*z/2)/math.sqrt(2*math.pi)
  scenarios=[('original_effective_N','original_rounded_P',raw,rawse,oldP),('original_effective_N','balanced_P_0_5',raw,rawse,.5),('total_N_derivative','original_rounded_P',new['tot'],new['tot_se'],oldP),('total_N_derivative','actual_original_case_fraction',new['tot'],new['tot_se'],actualP),('total_N_derivative','balanced_P_0_5',new['tot'],new['tot_se'],.5)]
  for arm,convention,h,s,Pvalue in scenarios:
   factor=K*K*(1-K)**2/(Pvalue*(1-Pvalue)*density*density);point=h*factor;ss=s*factor
   close(trait+'/'+arm+'/'+convention+'/presentation_h2_Z_invariance',point/ss,h/s,presentation=True)
   if arm=='original_effective_N' and convention=='original_rounded_P':close(trait+'/reported_h2_original',point,float(old['reported_h2_full_precision']),presentation=True);close(trait+'/reported_se_original',ss,float(old['reported_h2_se_full_precision']),presentation=True)
   presentations.append({'trait':trait,'N_arm':arm,'P_convention':convention,'historical_case_count':cases,'historical_control_count':controls,'total_N':cases+controls,'P_supplied':Pvalue,'K_historical_fixed':K,'raw_fitted_h2':h,'raw_fitted_SE':s,'liability_factor':factor,'presented_h2':point,'presented_SE':ss,'presented_CI95_lower':point-1.96*ss,'presented_CI95_upper':point+1.96*ss,'h2_Z_invariant':h/s,'intercept_invariant_under_presentation':float(old['intercept_full_precision']) if arm=='original_effective_N' else new['intercept'],'point_in_0_1':0<=point<=1,'historical_standalone_verdict':old['original_verdict'],'historical_exclusion_reason':old['original_qc_reason'],'presentation_only':True,'valid_population_h2_certified':False,'presentation_tolerance':1e-12,'no_new_QC_admission':True})
 out=T/'independent_sensitivity_62_full_precision_v4_6.tsv';pt=T/'independent_sensitivity_liability_presentations_v4_6.tsv';at=T/'independent_sensitivity_62_arithmetic_checks_v4_6.tsv';tsv(out,rows);tsv(pt,presentations);tsv(at,arithmetic)
 summary={}
 for group in ['lipid_two_step','binary_total_N']:
  selected=[r for r in rows if r['kind']=='rg' and r['job_id'].startswith(group)]
  summary[group]={'rows':len(selected),'baseline_intercept_both_pass':sum(r['baseline_pairwise_h2_intercept_le1_2_both'] for r in selected),'sensitivity_intercept_both_pass':sum(r['sensitivity_pairwise_h2_intercept_le1_2_both'] for r in selected),'intercept_boundaries_changed':sum(r['intercept_boundary_changed'] for r in selected),'baseline_h2_Z_both_pass':sum(r['baseline_pairwise_h2_Z_ge4_both'] for r in selected),'sensitivity_h2_Z_both_pass':sum(r['sensitivity_pairwise_h2_Z_ge4_both'] for r in selected),'h2_Z_boundaries_changed':sum(r['h2_Z_boundary_changed'] for r in selected),'max_abs_delta_rg':max(abs(r['delta_rg']) for r in selected),'max_abs_delta_rg_SE':max(abs(r['delta_rg_SE']) for r in selected),'warnings':sum(bool(r['sensitivity_warnings']) for r in selected)}
 receipt=R/'independent_sensitivity_whole62_receipt_v4_6.json'
 save(receipt,{'schema':'independent_actual_whole62_sensitivity_verification_v4_6','status':'ALL62_CAPTURE_DELETE_ARITHMETIC_AND_EXECUTION_BINDINGS_PASS_WITH_SCIENTIFIC_QUALIFICATIONS','checker_sha256':sha(__file__),'resource_plan_sha256':sha(R/'independent_sensitivity_whole62_resource_plan_v4_6.json'),'plan_sha256':PIN,'stage_sha256':STAGEPIN,'terminal_seal_sha256':SEALPIN,'checked_artifact_sha256':HASHES,'current_consumed_artifacts_verified':602,'baseline_receipt_count':len(stage['baseline_gate_receipt_sha256']),'stock_audit_commands':52,'paired_identity_instances':124,'fit_commands':26,'native_estimates':62,'all62_arithmetic_checks_pass':True,'stock_delete_arrays_exact':stock_arrays,'checks':CHECKS,'maximum_absolute_arithmetic_discrepancies':MAXERR,'descriptive_summary':summary,'standalone_liability_presentations':presentations,'table_sha256':{str(q):sha(q) for q in [out,pt,at]},'source_input_hash_maps_verified_from_native_receipts':True,'source_body_hashes_independently_reread':False,'baseline_stock_merge_is_historical_fit_instrumentation':False,'within_fit_relative_tolerance':1e-12,'within_fit_absolute_tolerance':1e-15,'within_fit_P_absolute_tolerance':1e-300,'presentation_relative_absolute_tolerance':1e-12,'no_tolerance_for_actual_sensitivity_differences':True,'historical396_q_and_exclusions_retained':True,'new_BH_families':0,'difference_P_values':0,'cross_fit_covariance_certified':False,'source_QC_resolved_by_sensitivity':False,'valid_population_h2_certified':False,'observed_RSS_bytes':usage(),'elapsed_seconds':time.monotonic()-START,'GWAS_body_reads':0,'reference_body_reads':0,'actual_fits':0,'actual_workers':0,'actual_mutex_acquisitions':0})
 print(json.dumps({'receipt_sha256':sha(receipt),'summary':summary,'checks':len(CHECKS),'RSS_bytes':usage(),'elapsed_seconds':time.monotonic()-START},indent=2))
if __name__=='__main__':main()
