#!/usr/bin/env python3
"""Automatic lexical screen of historical documents; does not claim peer review."""
import csv,hashlib,re
from pathlib import Path
from audit_evidence import ROOT,AREA,write_tsv
out=[]
for name in ['brain6','extensions/brain6','results/track_b']:
 for p in sorted((ROOT/name).rglob('*')):
  rel=str(p.relative_to(ROOT))
  if not p.is_file() or p.suffix not in ['.md','.txt','.json','.yml','.yaml'] or 'confirmatory_v2' in rel or 'work/' in rel:continue
  s=p.read_text(errors='replace');tags=[]
  for k,pat in [('PRIMARY_QC_FAIL',r'FAILED_QC_NOT_PROMOTED'),('PARTIAL_HISTORY',r'PARTIAL_FAMILY|19 regions|four.track'),('EXPLORATORY',r'EXPLORATORY|exploratory'),('REPLICATION_OR_CAUSAL_WORDING',r'replicat|causal|effector'),('ACCESS_HOLD',r'BLOCKED|NOT_RUN|unavailable')]:
   if re.search(pat,s,re.I):tags.append(k)
  out.append(dict(path=rel,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),semantic_tags=';'.join(tags),disposition='HISTORICAL_DOCUMENT_SCAN_NOT_A_CURRENT_CLAIM;authority=final_package_v1+continuation',method='AUTOMATED_LEXICAL_SCREEN_NOT_LINE_BY_LINE_HUMAN_REVIEW'))
write_tsv(AREA/'qc/document_claim_scan.tsv',out);print('Lexically screened historical documents/configs:',len(out))
