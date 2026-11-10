#!/usr/bin/env python3
import csv,hashlib,json,re
from pathlib import Path
import openpyxl
root=Path(__file__).resolve().parent.parent;p=root/'work/frailpol_participants.xlsx';gpath=root/'work/frailpol_gait_three_class.csv'
hashval=hashlib.md5(p.read_bytes()).hexdigest()
if hashval!='c438dc375281f2fb1e0afe08be89aed3':raise RuntimeError('Publisher checksum mismatch')
w=openpyxl.load_workbook(p,read_only=True,data_only=True);s=w.active;allrows=list(s.values);h=allrows[0];rows=[r for r in allrows[1:] if r[0] is not None and str(r[0]).isdigit()];d={str(r[0]):r for r in rows}
if len(d)!=len(rows):raise ValueError('duplicate participant IDs')
valid=[r for r in rows if all(r[i] in (0,1) for i in range(9,14)) and isinstance(r[14],(int,float))]
score=[r for r in rows if isinstance(r[14],(int,float)) and isinstance(r[15],(int,float))]
gait=list(csv.DictReader(gpath.open(newline='')));linked=[(a,d[a['participant_id']]) for a in gait if a['participant_id'] in d]
terms=re.compile(r'insomnia|sleep|bezsen|sennu|senność|zaburzenia snu',re.I)
result={'workbook_md5':hashval,'workbook_people':len(rows),'component_total_numeric':len(valid),'component_total_agree':sum(sum(r[9:14])==r[14] for r in valid),'target_total_numeric':len(score),'target_fried_mapping_agree':sum(r[15]==(0 if r[14]==0 else 1 if r[14] in (1,2) else 2) for r in score),'gait_rows':len(gait),'gait_id_match':len(linked),'gait_class_agree':sum(str(b[15])==a['Targets'] for a,b in linked),'gait_complete_components':sum(all(b[i] in (0,1) for i in range(9,14)) and isinstance(b[14],(int,float)) for a,b in linked),'gait_component_total_disagree':sum(all(b[i] in (0,1) for i in range(9,14)) and isinstance(b[14],(int,float)) and sum(b[9:14])!=b[14] for a,b in linked),'headers':[x for x in h if x is not None],'gait_derived_class_disagree':sum(all(b[i] in (0,1) for i in range(9,14)) and isinstance(b[14],(int,float)) and (0 if sum(b[9:14])==0 else 1 if sum(b[9:14])<3 else 2)!=b[15] for a,b in linked),'rows_with_explicit_sleep_terms':sum(any(isinstance(c,str) and terms.search(c) for c in r[1:]) for r in rows),'warning':'Counts from public dataset; exact terms not semantic review. No genotype or insomnia established. Five-component FI label is a physical-frailty count, not the multi-deficit frailty index; authors explain weakness is exhaustion.'}
(root/'outputs/stage6_frailpol_participant_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='headers'},indent=2))
