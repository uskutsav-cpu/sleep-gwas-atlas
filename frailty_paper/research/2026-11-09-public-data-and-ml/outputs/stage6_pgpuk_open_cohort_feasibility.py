#!/usr/bin/env python3
"""Aggregate-only PGP-UK API metadata feasibility scan; no individual data saved."""
import collections,hashlib,json,re,subprocess
from pathlib import Path
url='https://www.personalgenomes.org.uk/api/v1.3/all_participants'
raw=subprocess.run(['curl','-fsSL','--max-time','45',url],stdout=subprocess.PIPE,check=True).stdout
obj=json.loads(raw);assert isinstance(obj,dict)
terms={'positive_control_eye_colour':r'eye colou?r','sleep':r'\bsleep|\binsomnia|\bsleepless','frailty':r'\bfrail|\bfrailty','grip_strength':r'\bgrip|hand.?strength','walking':r'\bwalk|\bgait','fatigue':r'\bfatig|\btired|\bexhaust','weight_loss':r'weight (?:loss|lost)|los[st].*weight'}
counts={k:{'participants_with_matching_question':0,'participants_with_matching_answer_text':0} for k in terms};qentries=0;phen=0;variants=0;vcf=0
for participant in obj.values():
 forms=participant.get('phenotype') or [];phen+=bool(forms); files=participant.get('variant') or [];variants+=bool(files);vcf+=any('VCF' in str(x.get('submitted_format','')).upper() for x in files if isinstance(x,dict));found={k:[False,False] for k in terms}
 for form in forms:
  for item in form.get('data') or []:
   if not isinstance(item,dict):continue
   qentries+=1
   for key,pattern in terms.items():
    found[key][0]|=bool(re.search(pattern,str(item.get('question','')),re.I));found[key][1]|=bool(re.search(pattern,str(item.get('answer','')),re.I))
 for key,(question,answer) in found.items():
  counts[key]['participants_with_matching_question']+=question;counts[key]['participants_with_matching_answer_text']+=answer
result={'source':url,'source_sha256':hashlib.sha256(raw).hexdigest(),'participants_at_read':len(obj),'participants_with_nonempty_phenotype_metadata':phen,'participants_with_variant_metadata':variants,'participants_with_vcf_labelled_variant_metadata':vcf,'question_answer_entries_scanned':qentries,'case_insensitive_patterns':terms,'counts':counts,'scope':'Public API metadata only. No genomes or phenotypes linked or individually verified. Keyword hits are not diagnoses, non-hits are not proof of absence; participant overlap, file availability, variant coverage, ages and outcome data remain unverified. No individual identifiers or answers saved.'}
Path(__file__).with_name('stage6_pgpuk_open_cohort_feasibility.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['participants_at_read','participants_with_nonempty_phenotype_metadata','participants_with_variant_metadata','participants_with_vcf_labelled_variant_metadata','question_answer_entries_scanned','counts']},indent=2))
