#!/usr/bin/env python3
"""Write complete per-cell solver packets and preregister the bounded feasibility.

Freezing requires an external reviewer attestation; this script never calls a model.
"""
import argparse,hashlib,itertools,json
from datetime import datetime,timezone
from pathlib import Path

SEED=20261012
MODELS=['gpt-6.1-sol','gpt-6-sol','gpt-6-luna']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--public',type=Path,required=True);ap.add_argument('--private',type=Path,required=True);ap.add_argument('--freeze',action='store_true');ap.add_argument('--audit-attestation',type=Path);a=ap.parse_args()
 if (a.public/'FREEZE_RECEIPT.json').exists():raise ValueError('Already frozen; preserve original inventory and issue explicit versioned amendment')
 if a.freeze:
  if any((a.private/'model_responses').rglob('*.json')):raise ValueError('Model responses already exist: cannot declare pre-response freeze')
  if not a.audit_attestation or not a.audit_attestation.is_file():raise ValueError('Independent review attestation required before freeze')
  plan=json.loads((a.public/'execution_registration.json').read_text());plan['status']='FROZEN_BEFORE_FINAL_RESPONSES';plan['freeze_time_utc']=datetime.now(timezone.utc).isoformat();plan['independent_review_attestation_sha256']=sha(a.audit_attestation)
  # Verify all intended packets before changing registration state.
  for cell in plan['schedule']:
   p=a.private/cell['packet_relative_path']
   if sha(p)!=cell['packet_sha256']:raise ValueError('Packet hash drift: '+str(p))
  save(a.public/'execution_registration.json',plan);save(a.private/'execution_registration.json',plan)
  public_files=['PROTOCOL_AND_SAP.md','final_tasks.json','development_tasks.json','construction_manifest.json','execution_registration.json','build_benchmark.py','prepare_execution.py','score_benchmark.py','validate_benchmark.py','README.md']
  private_files=['GOLD_NOT_FOR_SOLVERS.json','condition_A_inputs.json','condition_B_inputs.json','condition_C_inputs.json']
  receipt={'freeze_time_utc':plan['freeze_time_utc'],'final_model_responses_at_freeze':0,'final_items':60,'development_items':6,'feasibility_final_item_ids':plan['selected_task_ids'],'planned_feasibility_cells':90,'main_60_execution':'PENDING','semantic_gold_certification':'Independent review is bounded; source validity and automatic gold semantic validity are not asserted','review_attestation_sha256':sha(a.audit_attestation),'public_files_sha256':{f:sha(a.public/f) for f in public_files},'private_files_sha256':{f:sha(a.private/f) for f in private_files},'packets_sha256':{f'packets/{c}/{item["task_id"]}.json':sha(a.private/f'packets/{c}/{item["task_id"]}.json') for c in 'ABC' for item in json.loads((a.private/f'condition_{c}_inputs.json').read_text())['items']}}
  save(a.public/'FREEZE_RECEIPT.json',receipt);save(a.private/'FREEZE_RECEIPT.json',receipt);print(json.dumps({'frozen':True,'time':receipt['freeze_time_utc'],'feasibility_cells':90}));return
 tasks=json.loads((a.public/'final_tasks.json').read_text());rank=lambda t:hashlib.sha256(f'{SEED}|feasibility|{t["task_id"]}'.encode()).hexdigest();strata=sorted({t['stratum'] for t in tasks});subset=[min((t for t in tasks if t['stratum']==s),key=rank) for s in strata];subset.sort(key=lambda t:hashlib.sha256(f'{SEED}|feasibility_order|{t["task_id"]}'.encode()).hexdigest())
 packets={}
 for c in 'ABC':
  data=json.loads((a.private/f'condition_{c}_inputs.json').read_text())
  for item in data['items']:
   p=a.private/'packets'/c/f'{item["task_id"]}.json';save(p,item);packets[c,item['task_id']]=p
 perms=list(itertools.permutations('ABC'));schedule=[]
 for pos,t in enumerate(subset):
  for mi,model in enumerate(MODELS):
   for ci,c in enumerate(perms[(pos*len(MODELS)+mi)%6]):
    p=packets[c,t['task_id']];schedule.append({'cell_index':len(schedule)+1,'task_position':pos+1,'condition_position':ci+1,'task_id':t['task_id'],'stratum':t['stratum'],'condition':c,'model_alias':model,'provider':'OpenAI','requested_reasoning_effort':'medium','requested_fork_turns':'none','packet_relative_path':f'packets/{c}/{t["task_id"]}.json','packet_sha256':sha(p),'packet_bytes':p.stat().st_size,'response_relative_path':f'model_responses/{model}/{c}/{t["task_id"]}.json','attempts_planned':1})
 plan={'status':'PREFREEZE_AWAITING_INDEPENDENT_REVIEW','seed':SEED,'experiment':'90-cell initial Codex-only feasibility; main60 pending','final_task_count':60,'selected_final_task_count':10,'selected_task_ids':[t['task_id'] for t in subset],'planned_cells':90,'models':[{'alias':m,'provider':'OpenAI','reasoning_effort_requested':'medium','exact_backend_version':'UNKNOWN','temperature':'UNKNOWN_UNCONTROLLABLE','sampling_seed':'UNKNOWN_UNCONTROLLABLE'} for m in MODELS],'technical_isolation':False,'instruction_only_isolation':True,'timing_and_tokens':'UNKNOWN unless host exposes','schedule':schedule,'complete_packet_max_bytes':max(p.stat().st_size for p in packets.values())}
 save(a.public/'execution_registration.json',plan);save(a.private/'execution_registration.json',plan)
 print(json.dumps({'subset':plan['selected_task_ids'],'planned_cells':90,'packets':len(packets),'max_bytes':plan['complete_packet_max_bytes']}))
if __name__=='__main__':main()
