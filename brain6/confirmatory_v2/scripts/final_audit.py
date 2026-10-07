#!/usr/bin/env python3
"""Read-only historical-file integrity and final continuation execution receipt."""
import csv,hashlib,json,subprocess
from pathlib import Path
from audit_evidence import ROOT,AREA
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
allowed={'tests/test_track_b_pleiotropy_results_v2.py','.github/workflows/ci.yml'};drift=[];missing=[];checked=0
for r in csv.DictReader((AREA/'qc/file_inventory.tsv').open(),delimiter='\t'):
 path=r.get('path') or r.get('file');expected=r.get('sha256')
 if not path or not expected:continue
 p=ROOT/path
 if not p.is_file():missing.append(path);continue
 checked+=1
 if sha(p)!=expected:drift.append(path)
unsafe=[p for p in drift if p not in allowed]
assert not unsafe,unsafe
r={'historical_file_hashes_checked':checked,'historical_missing_since_inventory':missing,'historical_hash_drift':drift,'permitted_engineering_changes':sorted(allowed),'scientific_output_drift':unsafe,'historical_scientific_outputs_preserved':not unsafe,'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'source_free_contracts':json.loads((AREA/'qc/contract_test_partition.json').read_text()),'manuscript_validation':json.loads((AREA/'qc/manuscript_validation.json').read_text()),'native_reproduction':'NOT_RUN_MISSING_ORIGINAL_INPUTS_REFERENCE_AND_PACKAGES','CI':'See qc/github_actions_status.json; local equivalent checks separately recorded','human_submission_approval':'NOT_PROVIDED'}
(AREA/'qc/final_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
files=[]
for p in sorted(AREA.rglob('*')):
 if p.is_file() and '__pycache__' not in str(p) and p.name not in ['DELIVERABLES.sha256','final_receipt.json'] and p.suffix not in ['.pyc']:
  files.append(sha(p)+'  '+str(p.relative_to(AREA)))
(AREA/'DELIVERABLES.sha256').write_text('\n'.join(files)+'\n')
print(json.dumps({k:r[k] for k in ['historical_file_hashes_checked','historical_hash_drift','historical_scientific_outputs_preserved','branch','head']}))
