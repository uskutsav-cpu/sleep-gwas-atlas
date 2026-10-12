#!/usr/bin/env python3
"""Score real structured responses against private condition-specific gold.

This checks frozen answer agreement, not independent semantic validity. Never
constructs model responses or fills unexecuted cells with synthetic answers.
"""
import argparse, hashlib, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

SEED=20261012

def exact(got,want,tolerance):
 if want is None: return got is None
 if isinstance(want,bool): return isinstance(got,bool) and got==want
 if isinstance(want,int): return isinstance(got,int) and not isinstance(got,bool) and got==want
 if isinstance(want,float):
  return isinstance(got,(int,float)) and not isinstance(got,bool) and math.isfinite(got) and abs(got-want)<=max(tolerance['absolute'],tolerance['relative']*abs(want))
 return isinstance(got,type(want)) and got==want

def structural(obj,schema):
 if not isinstance(obj,dict): return False
 if set(obj)!=set(schema['required']): return False
 for f,v in obj.items():
  prop=schema['properties'][f]; kinds=prop['type']; kinds=[kinds] if isinstance(kinds,str) else kinds
  if v is None:
   if 'null' not in kinds:return False
  elif isinstance(v,bool):
   if 'boolean' not in kinds:return False
  elif isinstance(v,str):
   if 'string' not in kinds:return False
  elif isinstance(v,int):
   if 'integer' not in kinds and 'number' not in kinds:return False
  elif isinstance(v,float):
   if 'number' not in kinds or not math.isfinite(v):return False
  else:return False
  if 'enum' in prop and v not in prop['enum']:return False
 return True

def unsafe(obj,gold):
 if not isinstance(obj,dict):return []
 risks=[]
 forbidden_true={'zero_sharing_established','independent_samples_supported','threshold_flag_alone_clears_inference','blank_signed_field_establishes_null','fully_independent_two_trait_replication','intervention_effect_established','calibrated_difference_test_available','biological_absence_established','global_novelty_established'}
 for f in forbidden_true:
  if f in gold and obj.get(f) is True:risks.append(f)
 if 'causal_direction' in gold and obj.get('causal_direction')=='ESTABLISHED':risks.append('causal_direction')
 if gold['status']=='WITHHELD_DIRECTION' and obj.get('signed_interpretation_withheld') is False:risks.append('direction_hold_cleared')
 if gold['status']=='WITHHELD_UNCERTAINTY' and obj.get('reliable_inference_withheld') is False:risks.append('uncertainty_hold_cleared')
 return sorted(risks)

def quantile(xs,p):
 xs=sorted(xs); at=p*(len(xs)-1); lo=int(at); hi=min(lo+1,len(xs)-1); return xs[lo]+(xs[hi]-xs[lo])*(at-lo)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--gold',type=Path,required=True);ap.add_argument('--tasks',type=Path,required=True);ap.add_argument('--responses',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--registration',type=Path,required=True);ap.add_argument('--bootstrap',type=int,default=10000);a=ap.parse_args()
 gold_data=json.loads(a.gold.read_text());gold={x['task_id']:x for x in gold_data['items'] if x['split']=='final'};tasks={x['task_id']:x for x in json.loads(a.tasks.read_text())}
 if set(tasks)!=set(gold):raise ValueError('Task/gold final inventory mismatch')
 registration=json.loads(a.registration.read_text()); selected=set(registration['selected_task_ids']); tasks={k:v for k,v in tasks.items() if k in selected}; gold={k:v for k,v in gold.items() if k in selected}; registered_models={m['alias'] for m in registration['models']}; planned={(c['model_alias'],c['condition'],c['task_id']) for c in registration['schedule']}; registered_cells={(c['model_alias'],c['condition'],c['task_id']):c for c in registration['schedule']}
 text=a.responses.read_text()
 try:raw=json.loads(text)
 except json.JSONDecodeError:raw=[json.loads(line) for line in text.splitlines() if line.strip()]
 if isinstance(raw,dict):raw=raw['responses']
 seen=set();scored=[];models=set();unevaluable=[]
 for r in raw:
  ident=(r['model_alias'],r['condition'],r['task_id'])
  if ident in seen:raise ValueError('Duplicate cell: '+repr(ident))
  seen.add(ident); model,c,tid=ident
  if c not in 'ABC' or tid not in gold:raise ValueError('Unknown condition/task')
  if r.get('provider')!='OpenAI':raise ValueError('Expected actual provider OpenAI, no relabelled providers')
  if ident not in planned:raise ValueError('Unregistered cell: '+repr(ident))
  if r.get('packet_sha256')!=registered_cells[ident]['packet_sha256']:raise ValueError('Missing or incorrect packet provenance: '+repr(ident))
  models.add(model)
  if r.get('execution_status','COMPLETED_RESPONSE')!='COMPLETED_RESPONSE':
   unevaluable.append({'model_alias':model,'condition':c,'task_id':tid,'execution_status':r.get('execution_status')});continue
  g=gold[tid];want=g['condition_answers'][c];obj=r.get('response');parse_ok=True
  if isinstance(obj,str):
   try:obj=json.loads(obj)
   except (json.JSONDecodeError,ValueError):parse_ok=False;obj=None
  schema_ok=parse_ok and structural(obj,tasks[tid]['response_schema']); field_correct={}
  for f,v in want.items():
   ok=isinstance(obj,dict) and f in obj and exact(obj[f],v,g['numeric_tolerance'])
   if not ok and c=='A' and isinstance(obj,dict) and f in obj:
    ok=any(exact(obj[f],alternative,g['numeric_tolerance']) for alternative in g.get('accepted_condition_A_alternatives',{}).get(f,[]))
   field_correct[f]=ok
  risk=unsafe(obj,want); numeric_fields=[f for f,v in want.items() if isinstance(v,(int,float)) and not isinstance(v,bool)]
  scored.append({'task_id':tid,'model_alias':model,'provider':'OpenAI','condition':c,'stratum':g['stratum'],'schema_valid':bool(schema_ok),'strict_correct':bool(schema_ok and all(field_correct.values())),'field_fraction':sum(field_correct.values())/len(field_correct),'field_correct':field_correct,'numeric_correct':sum(field_correct[f] for f in numeric_fields),'numeric_fields':len(numeric_fields),'unsafe_indicators':risk,'status_correct':field_correct['status'],'condition_A_metrics':None if c!='A' else {'exact_target_supported':want['status']!='INSUFFICIENT_TARGET_SOURCE_EVIDENCE','abstention_fields':[f for f,role in g['field_roles'].items() if role=='target_factual' and want[f] is None],'correct_abstention_fields':[f for f,role in g['field_roles'].items() if role=='target_factual' and want[f] is None and field_correct[f]],'unsupported_target_fields':[f for f,role in g['field_roles'].items() if role=='target_factual' and want[f] is None and isinstance(obj,dict) and f in obj and obj[f] is not None],'supported_conventional_fields':[f for f,role in g['field_roles'].items() if role=='conventional_factual' and want[f] is not None],'correct_conventional_fields':[f for f,role in g['field_roles'].items() if role=='conventional_factual' and want[f] is not None and field_correct[f]]},'pair_cluster':tasks[tid]['pair_cluster'],'source_clusters':tasks[tid]['source_clusters'],'response_sha256':hashlib.sha256(json.dumps(r.get('response'),sort_keys=True).encode()).hexdigest()})
 groups=defaultdict(list)
 for r in scored:groups[(r['model_alias'],r['condition'])].append(r)
 summaries=[]
 for (model,c),rows in sorted(groups.items()):
  n=len(rows);numeric_n=sum(r['numeric_fields'] for r in rows)
  summaries.append({'model_alias':model,'provider':'OpenAI','condition':c,'completed_final_cells':n,'expected_cells':len(selected),'strict_task_rate':sum(r['strict_correct'] for r in rows)/n,'mean_field_fraction':sum(r['field_fraction'] for r in rows)/n,'structural_validity_rate':sum(r['schema_valid'] for r in rows)/n,'status_accuracy':sum(r['status_correct'] for r in rows)/n,'applicable_numeric_accuracy':sum(r['numeric_correct'] for r in rows)/numeric_n if numeric_n else None,'unsafe_task_rate':sum(bool(r['unsafe_indicators']) for r in rows)/n,'strata':{s:{'n':len(rs),'strict_task_rate':sum(x['strict_correct'] for x in rs)/len(rs)} for s,rs in ((s,[r for r in rows if r['stratum']==s]) for s in sorted({r['stratum'] for r in rows}))}})
 index={(r['model_alias'],r['condition'],r['task_id']):r for r in scored}; contrasts=[]
 for model in sorted(models):
  paired=[tid for tid in sorted(tasks) if (model,'B',tid) in index and (model,'C',tid) in index]
  if not paired:continue
  differences={tid:float(index[model,'C',tid]['strict_correct'])-float(index[model,'B',tid]['strict_correct']) for tid in paired}
  contrasts.append({'model_alias':model,'paired_tasks':len(paired),'strict_task_rate_difference_C_minus_B':sum(differences.values())/len(differences),'field_fraction_difference_C_minus_B':sum(index[model,'C',tid]['field_fraction']-index[model,'B',tid]['field_fraction'] for tid in paired)/len(paired),'task_differences':differences})
 # Only taskwise complete B/C observations are pooled. Missingness is disclosed,
 # never silently imputed and never treated as an incorrect scientific response.
 differences={tid:sum(x['task_differences'][tid] for x in contrasts if tid in x['task_differences'])/sum(tid in x['task_differences'] for x in contrasts) for tid in tasks if any(tid in x['task_differences'] for x in contrasts)}
 dependence=None
 if differences and len(selected)==60 and len(differences)==60:
  sleeps=sorted({tasks[t]['source_clusters'][0] for t in differences});outcomes=sorted({tasks[t]['source_clusters'][1] for t in differences});rng=random.Random(SEED);samples=[]
  for _ in range(a.bootstrap):
   sw={s:rng.expovariate(1) for s in sleeps};ow={o:rng.expovariate(1) for o in outcomes};weights={t:sw[tasks[t]['source_clusters'][0]]*ow[tasks[t]['source_clusters'][1]] for t in differences};samples.append(sum(weights[t]*differences[t] for t in differences)/sum(weights.values()))
  leave=[]
  for side,clusters in [(0,sleeps),(1,outcomes)]:
   for cluster in clusters:
    vals=[v for t,v in differences.items() if tasks[t]['source_clusters'][side]!=cluster]
    if vals:leave.append(sum(vals)/len(vals))
  dependence={'taskwise_pooled_C_minus_B':sum(differences.values())/len(differences),'task_count':len(differences),'fixed_executed_model_count':len(contrasts),'sleep_source_clusters':len(sleeps),'outcome_source_clusters':len(outcomes),'method':'two-way source-weighted Exp(1) bootstrap sensitivity; not a provider-population CI','replicates':a.bootstrap,'seed':SEED,'percentile_95_interval':[quantile(samples,.025),quantile(samples,.975)] if samples else None,'leave_one_source_out_range':[min(leave),max(leave)] if leave else None}
 a_rows=[r['condition_A_metrics'] for r in scored if r['condition']=='A']; abstention_n=sum(len(x['abstention_fields']) for x in a_rows);conv_n=sum(len(x['supported_conventional_fields']) for x in a_rows)
 a_metrics={'evaluated_A_cells':len(a_rows),'exact_target_source_coverage':sum(x['exact_target_supported'] for x in a_rows)/len(a_rows) if a_rows else None,'required_target_abstention_fields':abstention_n,'valid_abstention_field_rate':sum(len(x['correct_abstention_fields']) for x in a_rows)/abstention_n if abstention_n else None,'supported_conventional_factual_fields':conv_n,'supported_conventional_factual_yield':sum(len(x['correct_conventional_fields']) for x in a_rows)/conv_n if conv_n else None,'unsupported_target_substitution_task_rate':sum(bool(x['unsupported_target_fields']) for x in a_rows)/len(a_rows) if a_rows else None}
 report={'A_answerability_and_abstention':a_metrics,'registered_models':sorted(registered_models),'planned_cells':len(planned),'observed_cell_records':len(seen),'unobserved_planned_cells':len(planned-seen),'unevaluable_cells':unevaluable,'subset_inference_policy':'No confidence interval for ten-item feasibility; no general interface/provider superiority claim.','scientific_gold_status':sorted({g['semantic_gold_status'] for g in gold.values()}),'score_scope':'Mechanical agreement with frozen condition-specific structured gold; semantic validity requires independent scientific adjudication. Exploratory Codex-only agent workflow.','completed_cells':len(scored),'planned_cells_for_observed_models':3*len(selected)*len(models),'provider':'OpenAI','model_aliases':sorted(models),'cross_provider_generalization':False,'summaries':summaries,'paired_B_C':contrasts,'dependence_sensitivity':dependence,'A_policy':'Do not compare pooled A/B/C accuracy as representation effect. Report supported factual yield, source coverage and correct abstention separately.','source_hashes':{'gold':hashlib.sha256(a.gold.read_bytes()).hexdigest(),'tasks':hashlib.sha256(a.tasks.read_bytes()).hexdigest(),'responses':hashlib.sha256(a.responses.read_bytes()).hexdigest()}}
 a.out.mkdir(parents=True,exist_ok=True);(a.out/'scored_cells.json').write_text(json.dumps(scored,indent=2)+'\n');(a.out/'score_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'completed_cells':len(scored),'out':str(a.out)}))
if __name__=='__main__':main()
