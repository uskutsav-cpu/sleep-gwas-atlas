#!/usr/bin/env python3
"""Independently reconstruct factual parity and check draft integrity, not validity."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_schema(v,s):
 t=s['type']
 if t=='object':
  assert isinstance(v,dict) and set(v)==set(s['required'])
  for k,x in v.items():check_schema(x,s['properties'][k])
 elif t=='array':
  assert isinstance(v,list)
  if s.get('uniqueItems'):assert len(v)==len({json.dumps(x,sort_keys=True) for x in v})
  for x in v:check_schema(x,s['items'])
 elif t=='string':assert isinstance(v,str) and ('enum' not in s or v in s['enum'])
 elif t=='number':assert isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
 elif t=='boolean':assert isinstance(v,bool)
 else:raise AssertionError(t)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--public',type=Path,required=True);ap.add_argument('--private',type=Path,required=True);ap.add_argument('--repo',type=Path,required=True);a=ap.parse_args();tasks=json.loads((a.public/'DRAFT_TASKS.json').read_text());g=json.loads((a.private/'DRAFT_GOLD_NOT_FOR_SOLVERS.json').read_text())['items'];m=json.loads((a.public/'DRAFT_CONSTRUCTION_MANIFEST.json').read_text());assert len(tasks)==len(g)==6
 master={c:{x['task_id']:x for x in json.loads((a.private/f'DRAFT_condition_{c}_inputs.json').read_text())['items']} for c in ['B','C']};prior=a.repo/'sleep_interpretation_validation_20261012/benchmark';old=json.loads((prior/'final_tasks.json').read_text())+json.loads((prior/'development_tasks.json').read_text());excluded={tuple(x['record_key']) for x in old};sourcekeys=set();packet_total=0
 for task,gold in zip(tasks,g):
  tid=task['task_id'];assert tid==gold['task_id'];assert gold['condition_answers']['B']==gold['condition_answers']['C'];check_schema(gold['condition_answers']['B'],task['response_schema']);b=master['B'][tid]['evidence'];c=master['C'][tid]['evidence'];assert b['common']==c['common'];bt=b['separate_factual_tables'];cards=c['source_cards'];tables={'estimates':[x['estimate'] for x in cards]}
  for table,cardfields,keyfields in [('sleep_constructs',['sleep_construct'],['trait_id']),('outcome_phenotypes',['outcome_phenotype'],['analysis_family','trait_id']),('heritability',['sleep_heritability','outcome_heritability'],['analysis_family','trait_id'])]:
   unique={tuple(x[f][k] for k in keyfields):x[f] for x in cards for f in cardfields if f in x};tables[table]=[unique[tuple(x[k] for k in keyfields)] for x in bt[table]]
  tables.update(c.get('other_factual_tables',{}))
  if 'conventional_source_cards' in c:
   hc=c['conventional_source_cards'];hm={x[f]['id']:x[f] for x in hc for f in ['id1_metadata','id2_metadata']};tables['hga_gc']=[x['gc'] for x in hc];tables['hga_metadata']=[hm[i] for i in sorted(hm,key=int)]
  assert tables==bt
  for row in bt['estimates']:
   rk=tuple(row[f] for f in ['analysis_family','sleep_trait','outcome_trait']);assert rk not in excluded;sourcekeys.add(rk)
  for condition in ['B','C']:
   p=a.private/'packets'/condition/f'{tid}.json';assert json.loads(p.read_text())==master[condition][tid];assert sha(p)==m['private_packets_sha256'][f'{condition}/{tid}'];packet_total+=p.stat().st_size
 # Recompute each original-file pointer independently. HGA compressed pointers
 # are reviewed from original cache separately; do not copy its statistics here.
 release=Path('/Volumes/Extreme SSD/sleep_interpretation_validation_20261012/source_qualified_candidate');checked=0
 for gold in g:
  for pointer in gold['row_field_provenance']:
   if pointer['file'].endswith('.gz'):continue
   p=(a.repo/pointer['file']) if pointer['file'].startswith('config/') else release/pointer['file'];assert sha(p)==pointer['sha256']
   with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
   r=rows[pointer['line']-2]
   if 'row_key' in pointer:assert all(r[k]==v for k,v in pointer['row_key'].items())
   checked+=1
 assert len(sourcekeys)==22
 frozen=json.loads((prior/'FREEZE_RECEIPT.json').read_text())
 for f,h in frozen['public_files_sha256'].items():assert sha(prior/f)==h
 report={'status':'DRAFT_MECHANICAL_CHECKS_ONLY','draft_items':6,'matched_condition_packets':12,'full_factual_object_parity':'6/6','gold_closed_schema':'6/6','new_atlas_pair_records':22,'exact_prior66_pair_overlap':0,'checked_original_TSV_pointers':checked,'original60_freeze_unchanged':True,'serialized_packet_bytes_total':packet_total,'max_packet_bytes':m['max_packet_bytes'],'model_responses':0,'human_scientific_certification':'PENDING','semantic_gold_validity':'NOT_CERTIFIED_BY_CODE'};print(json.dumps(report,indent=2));(a.public/'DRAFT_MECHANICAL_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n');(a.private/'DRAFT_MECHANICAL_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
