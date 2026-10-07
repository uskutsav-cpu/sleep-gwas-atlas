#!/usr/bin/env python3
"""Fail-closed manuscript claim and reference checks; no native inference claim."""
import csv,hashlib,json,re
from pathlib import Path
from audit_evidence import ROOT,AREA

def validate_text(text):
 issues=[]
 # These phrases are unsupported regardless of numerical output or prose placement.
 for forbidden in ['we identified a novel shared causal variant','we established independent two-trait replication','canonical LAVA passed QC','all LD matrices passed','a causal effector gene was confirmed']:
  if forbidden.lower() in text.lower():issues.append('UNSUPPORTED_POSITIVE_CLAIM:'+forbidden)
 for required in ['FAILED\\_QC\\_NOT\\_PROMOTED','3,720','873','8,465','1.1037','single-causal-variant','No assessed LD matrix','Independent two-trait locus replication was not established','full native-input reproduction remains incomplete']:
  if required.lower() not in text.lower():issues.append('MISSING_REQUIRED_LIMIT_OR_COUNT:'+required)
 return issues

def main():
 p=AREA/'MANUSCRIPT/brain6_manuscript.tex';text=p.read_text();issues=validate_text(text)
 claim=list(csv.DictReader((AREA/'CLAIM_TO_EVIDENCE.tsv').open(),delimiter='\t'))
 for c in claim:
  if hashlib.sha256((ROOT/c['exact_input']).read_bytes()).hexdigest()!=c['input_sha256']:issues.append('CLAIM_SOURCE_DRIFT:'+c['claim_id'])
 refs=list(csv.DictReader((AREA/'MANUSCRIPT/REFERENCES_VERIFIED.tsv').open(),delimiter='\t'));keys={r['key'] for r in refs};used=set()
 for group in re.findall(r'\\cite\{([^}]+)\}',text):used.update(group.split(','))
 if used-keys:issues.append('UNVERIFIED_CITATIONS:'+','.join(used-keys))
 a=text.split('\\section*{Abstract}')[1].split('\\textbf{Keywords:}')[0];maintext=text.split('\\section*{Background}')[1].split('\\section*{Abbreviations}')[0]
 def words(s):
  s=re.sub(r'\\(?:cite|href|texttt)\{[^}]*\}','',s);s=re.sub(r'\\[A-Za-z]+\*?','',s);return len(re.findall(r"\b[\w]+(?:[’'-][\w]+)*\b",s))
 n=words(a);mn=words(maintext)
 if n>350:issues.append('GM_ABSTRACT_TOO_LONG')
 if mn>3500:issues.append('CONSERVATIVE_MP_MAIN_TEXT_LIMIT_EXCEEDED')
 report={'status':'PASS' if not issues else 'FAIL','issues':issues,'structured_abstract_words_approx':n,'main_words_approx_including_table':mn,'verified_reference_count':len(refs),'cited_reference_count':len(used),'claims_checked':len(claim),'scope':'Exact allowed claims, limits, source hashes and indexed reference IDs; does not verify native analyses, author approval or scientific merit'}
 (AREA/'qc/manuscript_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return 0 if not issues else 1
if __name__=='__main__':raise SystemExit(main())
