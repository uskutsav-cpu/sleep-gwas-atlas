#!/usr/bin/env python3
"""Mechanical inventory/parity/provenance checks; not semantic certification."""
import argparse,hashlib,importlib.util,json
from collections import Counter
from pathlib import Path

def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--public',type=Path,required=True);ap.add_argument('--private',type=Path,required=True);a=ap.parse_args()
 build=module(a.public/'build_benchmark.py','benchmark_builder');score=module(a.public/'score_benchmark.py','benchmark_scorer')
 final=json.loads((a.public/'final_tasks.json').read_text());dev=json.loads((a.public/'development_tasks.json').read_text());tasks={t['task_id']:t for t in final+dev};gold=json.loads((a.private/'GOLD_NOT_FOR_SOLVERS.json').read_text())['items'];g={x['task_id']:x for x in gold};manifest=json.loads((a.public/'construction_manifest.json').read_text());registration=json.loads((a.public/'execution_registration.json').read_text())
 assert len(final)==60 and len(dev)==6 and len(tasks)==66 and set(tasks)==set(g)
 assert set(Counter(t['stratum'] for t in final).values())=={6} and len({t['stratum'] for t in final})==10
 keys=[tuple(t['record_key']) for t in final+dev];assert len(keys)==len(set(keys));assert not set(keys)&build.EXCLUDED
 packets={c:{i['task_id']:i for i in json.loads((a.private/f'condition_{c}_inputs.json').read_text())['items']} for c in 'ABC'}
 maxchars=maxbytes=0
 for tid,t in tasks.items():
  assert build.rejoin(packets['B'][tid]['evidence'])==packets['C'][tid]['evidence']
  assert g[tid]['condition_answers']['B']==g[tid]['condition_answers']['C']
  assert set(g[tid]['field_roles'])==set(g[tid]['condition_answers']['A'])
  for c in 'ABC':
   p=packets[c][tid];assert p['user_prompt']==t['prompt'] and p['response_schema']==t['response_schema'];assert score.structural(g[tid]['condition_answers'][c],t['response_schema']);assert not score.unsafe(g[tid]['condition_answers'][c],g[tid]['condition_answers'][c]);assert 'condition_answers' not in p and 'field_provenance' not in p
   ip=a.private/'packets'/c/f'{tid}.json';assert json.loads(ip.read_text())==p;maxchars=max(maxchars,len(json.dumps(p)));maxbytes=max(maxbytes,ip.stat().st_size)
 assert score.exact(None,None,{'absolute':1e-6,'relative':1e-5});assert not score.exact(0,None,{'absolute':1e-6,'relative':1e-5});assert not score.exact(True,1,{'absolute':1e-6,'relative':1e-5})
 assert not score.unsafe({'status':'AVAILABLE','reliable_inference_withheld':True},{'status':'WITHHELD_UNCERTAINTY','reliable_inference_withheld':True})
 assert score.unsafe({'status':'WITHHELD_UNCERTAINTY','reliable_inference_withheld':False},{'status':'WITHHELD_UNCERTAINTY','reliable_inference_withheld':True})==['uncertainty_hold_cleared']
 assert len(registration['schedule'])==90 and len({(x['task_id'],x['condition'],x['model_alias']) for x in registration['schedule']})==90
 selected={min([t for t in final if t['stratum']==s],key=lambda t:hashlib.sha256(f'20261012|feasibility|{t["task_id"]}'.encode()).hexdigest())['task_id'] for s in {t['stratum'] for t in final}}
 assert set(registration['selected_task_ids'])==selected
 for cell in registration['schedule']:assert sha(a.private/cell['packet_relative_path'])==cell['packet_sha256']
 for f,h in manifest['public_artifact_sha256'].items():assert sha(a.public/f)==h
 for f,h in manifest['private_artifact_sha256'].items():assert sha(a.private/f)==h
 if (a.public/'FREEZE_RECEIPT.json').exists():
  freeze=json.loads((a.public/'FREEZE_RECEIPT.json').read_text())
  for f,h in freeze['public_files_sha256'].items():assert sha(a.public/f)==h
  for f,h in freeze['private_files_sha256'].items():assert sha(a.private/f)==h
 report={'check':'mechanical_only_not_semantic_certification','final_items':60,'development_items':6,'unique_records':66,'B_C_full_object_parity':'66/66','gold_schema_validation':'198/198','planned_feasibility_cells':90,'feasibility_task_ids':registration['selected_task_ids'],'max_packet_bytes':maxbytes,'max_packet_approx_json_characters':maxchars,'A_final_exact_target_support':sum(x['condition_answers']['A']['status']!='INSUFFICIENT_TARGET_SOURCE_EVIDENCE' for x in gold if x['split']=='final'),'A_final_supported_conventional_fields':sum(sum(x['field_roles'][f]=='conventional_factual' and v is not None for f,v in x['condition_answers']['A'].items()) for x in gold if x['split']=='final'),'semantic_gold_status':'provisional; independent human scientific review pending'}
 print(json.dumps(report,indent=2));(a.private/'mechanical_validation_report.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
