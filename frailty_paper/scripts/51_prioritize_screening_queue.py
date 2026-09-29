#!/usr/bin/env python3
"""Create a transparent, non-exclusion keyword-priority list for human screening."""
from __future__ import annotations
import argparse, csv, hashlib, json, re, sys
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('frailty_paper/review/screening/title_abstract_queue.tsv')
PROTOCOL=Path('frailty_paper/review/protocol.md')
DEFAULT_OUT=Path('frailty_paper/review/screening/title_abstract_triage_priority.tsv')
RULES={
 'sleep_exposure': (3, [r'\bsleep\b',r'insomnia',r'chronotyp\w*',r'circadian',r'nap(?:ping)?',r'sleepiness',r'snoring',r'apn[oe]a',r'\bosa\b',r'sleep.?disordered breathing',r'sleep duration',r'sleep quality',r'sleep efficiency',r'sleep fragmentation']),
 'primary_frailty': (4, [r'frailt\w*',r'\bFried\b',r'clinical frailty scale',r'frailty index',r'frailty phenotype',r'pre.?frail']),
 'secondary_aging_outcome': (1, [r'sarcopen\w*',r'grip strength',r'walking speed',r'gait speed',r'mobility',r'physical activit\w*',r'disabilit\w*',r'fall\w*',r'healthspan',r'longevity',r'cognit\w*',r'mortality']),
 'observational_design': (1, [r'cohort',r'cross.?sectional',r'case.?control',r'longitudinal',r'population.?based',r'prospective',r'retrospective']),
 'genetic_method': (1, [r'\bGWAS\b',r'genetic correlation',r'Mendelian randomi[sz]ation',r'\bMR\b',r'genetic association',r'fine.?mapping',r'colocali[sz]ation',r'QTL'])
}
TITLE_BONUS=1

def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def compile_rules():return {k:(weight,[re.compile(x,re.I) for x in pats]) for k,(weight,pats) in RULES.items()}
def hits(text:str,rules)->dict[str,list[str]]:
 out={}
 for group,(_,patterns) in rules.items():
  found=[]
  for pat in patterns:
   m=pat.search(text)
   if m:
    term=m.group(0).lower()
    if term not in found:found.append(term)
  if found:out[group]=found
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo',type=Path,default=ROOT)
 ap.add_argument('--output',type=Path,default=DEFAULT_OUT)
 args=ap.parse_args();repo=args.repo.resolve();src=repo/SOURCE;protocol=repo/PROTOCOL;out=(repo/args.output if not args.output.is_absolute() else args.output).resolve()
 if not src.is_file() or not protocol.is_file():raise SystemExit('Queue or frozen protocol is missing')
 rules=compile_rules(); rows=[]
 with src.open(newline='',encoding='utf-8') as f:
  reader=csv.DictReader(f,delimiter='\t')
  required={'screening_id','pmid','doi','year','title','abstract','query_id','duplicate_count'}
  if not reader.fieldnames or not required.issubset(reader.fieldnames):raise SystemExit('Screening queue schema mismatch')
  for row in reader:
   title=row.get('title','') or ''; abstract=row.get('abstract','') or ''; joined=title+'\n'+abstract
   groups=hits(joined,rules); score=sum(rules[g][0] for g in groups)
   title_groups=hits(title,rules)
   title_bonus=sum(TITLE_BONUS for g in title_groups if g in ('sleep_exposure','primary_frailty'))
   score+=title_bonus
   direct=bool('sleep_exposure' in groups and 'primary_frailty' in groups)
   rows.append({'screening_id':row['screening_id'],'pmid':row.get('pmid',''),'doi':row.get('doi',''),'year':row.get('year',''),'title':title,'query_id':row.get('query_id',''),'triage_score':score,'direct_sleep_frailty_term_overlap':'YES' if direct else 'NO','matched_term_groups':';'.join(groups),'matched_terms':';'.join(f'{g}={"|".join(v)}' for g,v in groups.items()),'abstract_missing':'YES' if not abstract.strip() else 'NO','title_or_abstract_only':'YES'})
 if len({r['screening_id'] for r in rows})!=len(rows):raise SystemExit('Duplicate screening_id in source queue; refusing triage build')
 rows.sort(key=lambda r:(-r['triage_score'],r['screening_id']))
 fields=['priority_rank','screening_id','pmid','doi','year','title','query_id','triage_score','direct_sleep_frailty_term_overlap','matched_term_groups','matched_terms','abstract_missing','title_or_abstract_only']
 out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for rank,row in enumerate(rows,1):writer.writerow({'priority_rank':rank,**row})
 counts={}
 for row in rows:
  label=';'.join(g for g in RULES if g in row['matched_term_groups'].split(';')) or 'no_keyword_match'
  counts[label]=counts.get(label,0)+1
 manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'script':'frailty_paper/scripts/51_prioritize_screening_queue.py','script_sha256':sha(Path(__file__).resolve()),'command':'python frailty_paper/scripts/51_prioritize_screening_queue.py'+(' '+' '.join(sys.argv[1:]) if len(sys.argv)>1 else ''),'source_queue':SOURCE.as_posix(),'source_queue_sha256':sha(src),'frozen_protocol':PROTOCOL.as_posix(),'protocol_sha256':sha(protocol),'output':out.relative_to(repo).as_posix() if out.is_relative_to(repo) else out.as_posix(),'output_sha256':sha(out),'records':len(rows),'direct_sleep_frailty_term_overlap':sum(r['direct_sleep_frailty_term_overlap']=='YES' for r in rows),'missing_abstracts':sum(r['abstract_missing']=='YES' for r in rows),'priority_rule':{'method':'deterministic keyword grouping; binary group weights; title hit adds one point for sleep and/or primary frailty group','weights':{k:v[0] for k,v in RULES.items()},'groups':{k:v[1] for k,v in RULES.items()},'decisions':'NO screening decisions made; all records retained in priority list and source queue unchanged','purpose':'ordering assistance only; reviewers must independently screen every record against protocol','limitations':['Keyword matches are not eligibility judgments.','No record is excluded, included, or marked unclear by this script.','No age, population, design, or phenotype eligibility is inferred from missing terms.','Abstract-missing records remain in the list and are flagged.','This ranked list may bias reviewer attention; reviewers should preserve independent decisions and can ignore the order.']},'matched_group_counts':counts}
 manifest_path=out.with_suffix('.provenance.json');manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8')
 print(json.dumps({'records':len(rows),'direct_sleep_frailty_term_overlap':manifest['direct_sleep_frailty_term_overlap'],'missing_abstracts':manifest['missing_abstracts'],'output':str(out),'provenance':str(manifest_path),'queue_unchanged':True}))
if __name__=='__main__':main()
