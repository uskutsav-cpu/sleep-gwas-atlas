#!/usr/bin/env python3
"""Empirical source/accounting checks, not synthetic biology or result selection."""
import csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];AREA=ROOT/'brain6/translational_psychiatry_research_v1/research_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def main():
 checks=[]
 def check(name,ok,details):checks.append(dict(check=name,status='PASS' if ok else 'FAIL',details=details))
 master=read(AREA/'NOVELTY_MASTER.tsv');p=[r for r in master if r['candidate_locus_id'].startswith('placo_')]
 check('Complete immutable protected candidate accounting',len(p)==25 and len({r['candidate_locus_id'] for r in p})==25,'25 pair-specific IDs; no merged display interval substitution')
 check('Observed primary prior-art count',sum(r['classification']=='KNOWN_SAME_PAIR_REGION' for r in p)==22,'22 of25; remaining3 longsleep__parkinson unresolved, not novel')
 check('No absence-based novelty label',all('NOVEL' not in r['classification'].replace('NO_NOVELTY_INFERENCE','') for r in master),'No unresolved row promoted to novel')
 queries=read(AREA/'novelty/query_completion.tsv')
 check('All5 exact indexed queries fully paginated',len(queries)==5 and all(r['complete_for_query']=='True' and r['source_hit_count']==r['unique_records_retrieved'] for r in queries),'Complete exact queries only; not global exhaustive novelty proof')
 sources=read(AREA/'novelty/source_build_checks.tsv')
 check('All queried lead mappings coordinate consistent',all(x['coordinate_match']=='True' for x in sources),str(len(sources))+' selected official Ensembl37 mapping checks; absent mappings remain unverified')
 signals=read(AREA/'novelty/published_candidate_leads.tsv')
 for study,path in [('Lin_2026','lin_supp1.xlsx'),('Jia_2025','jia_supp1.docx')]:
  file=ROOT/'work/tp-literature-20261007'/path;check(study+' source supplement hashes',all(r['source_sha256']==sha(file) for r in signals if r['study']==study),sha(file))
 frozen=AREA/'sources/FIXED_VARIANT_FINGEN_MANIFEST.json';m=json.loads(frozen.read_text())
 check('Immutable freeze code/protocol hashes',all(sha(AREA/'sources'/n)==m['source_and_code_sha256'][n] for n in ['fixed_variant_finngen.py','FIXED_VARIANT_FINGEN_PROTOCOL.md']),sha(frozen))
 rec=json.loads((AREA/'sources/fixed_variant_acquisition_receipts.json').read_text())
 check('All3 selected raw responses acquired after freeze',len(rec)==3 and all(r['retrieved_utc']>=m['frozen_utc'] and r['status']=='RETRIEVED' and sha(ROOT/r['raw_path'])==r['sha256'] for r in rec),'Three exact point objects; hashes and dates verified')
 slots=read(AREA/'sources/fixed_variant_results_corrected.tsv');native=read(AREA/'sources/fixed_variant_native_fields_corrected.tsv')
 check('All4 planned slots retain fixed multiplicity',len(slots)==4 and {r['slot_id'] for r in slots}=={r['slot_id'] for r in m['slots']} and all(float(r['alpha_per_slot'])==.0125 for r in slots),'Two estimated threshold successes, one estimated nonpass, one NOT_ESTIMATED; no shrinkage')
 missing=[r for r in slots if r['decision']=='NOT_ESTIMATED']
 check('Missing P not imputed',len(missing)==1 and missing[0]['slot_id']=='A_INS_ADHD' and missing[0]['conjunction_p']==missing[0]['bonferroni_p']=='','A ADHD null statistics retained; no biological null')
 check('Seven planned native point rows retained',len(native)==7 and sum(x['status']=='ESTIMATED' for x in native)==6,'A3/B2/C2 endpoints; one missing placeholder')
 rep=read(AREA/'REPLICATION_MASTER.tsv');check('No unverified independent two-trait replication',all(r['classification']!='INDEPENDENT_TWO_TRAIT' for r in rep),'Participant-level/source-control linkage not available')
 graph=read(AREA/'sources/ADHD_DISCOVERY_COHORTS.tsv');check('Independent source manifests total correctly',len(graph)==13 and sum(int(x['cases']) for x in graph)==38691 and sum(int(x['controls']) for x in graph)==186843 and all(x['cross_source_count_match']=='True' for x in graph),'13 cohort names,26 count matches in two primary manifests')
 with (AREA/'novelty/source_accounting_checks.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
 print(json.dumps(dict(checks=len(checks),passed=sum(x['status']=='PASS' for x in checks),failed=sum(x['status']=='FAIL' for x in checks)),indent=2))
 if any(x['status']=='FAIL' for x in checks):sys.exit(1)
if __name__=='__main__':main()
