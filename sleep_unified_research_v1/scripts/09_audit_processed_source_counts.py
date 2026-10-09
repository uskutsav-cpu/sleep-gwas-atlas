#!/usr/bin/env python3
"""Low-memory direct count audit of all recovered extension input templates."""
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import time

PACKAGE=Path(__file__).resolve().parents[1]
DATA=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
def main():
    with (PACKAGE/'tables/archived_extension_recovery.tsv').open() as f:
        admitted={r['path']:r for r in csv.DictReader(f,delimiter='\t')}
    out=PACKAGE/'tables/extension_processed_source_counts_v1.tsv'
    if out.exists():raise SystemExit('PRIOR_COUNT_AUDIT_PRESERVED')
    rows=[];started=time.time()
    paths=sorted(p for p in (DATA/'data/munged').glob('*.sumstats.gz') if not p.name.startswith('._'))
    if len(paths)!=100:raise SystemExit('LOCKED100_SOURCE_COUNT_FAILED')
    with out.open('x',newline='') as f:
        fields=['trait_id','source_path','source_sha256','expected_sha256','total_template_rows','finite_N_Z_rows','missing_or_nonfinite_N_Z_rows','nonpositive_N_rows','h2_log_input_count','h2_log_count_matches_finite_rows','historical_rg_template_count','rg_count_matches_total_rows']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for i,path in enumerate(paths,1):
            claim=admitted[str(path)]
            actual=sha(path)
            if claim['status']!='EXACT_RECEIPT_HASH_RECOVERED' or actual!=claim['expected_sha256']:raise SystemExit('HASH_GATE_FAILED')
            n=finite=nonpositive=0
            with gzip.open(path,'rt') as source:
                header=source.readline().rstrip().split('\t')
                if header!=['SNP','A1','A2','Z','N']:raise SystemExit('MUNGED_SCHEMA_CHANGED')
                for line in source:
                    parts=line.rstrip().split('\t');n+=1
                    try:
                        z=float(parts[3]);sample=float(parts[4])
                        good=math.isfinite(z) and math.isfinite(sample)
                    except (ValueError,IndexError):good=False
                    finite+=good
                    if good and sample<=0:nonpositive+=1
            trait=path.name.replace('.sumstats.gz','')
            log=DATA/'logs/h2'/f'h2_{trait}.log'
            h2n=int(re.search(r'Read summary statistics for (\d+) SNPs',log.read_text())[1])
            # All twelve retained rg logs report the same complete HM3 template count.
            rgcounts=[]
            for rglog in sorted(p for p in (DATA/'logs/rg').glob('*.log') if not p.name.startswith('._')):
                text=rglog.read_text()
                match=re.search(r'Reading summary statistics from [^\n]*'+re.escape(trait)+r'\.sumstats\.gz \.\.\.\s*\nRead summary statistics for (\d+) SNPs',text)
                if match is None:raise SystemExit('RG_INPUT_COUNT_MISSING')
                rgcounts.append(int(match[1]))
            if len(rgcounts)!=12 or len(set(rgcounts))!=1:raise SystemExit('RG_TEMPLATE_COUNT_NOT_COMMON')
            row=dict(trait_id=trait,source_path=str(path),source_sha256=actual,expected_sha256=claim['expected_sha256'],total_template_rows=n,finite_N_Z_rows=finite,missing_or_nonfinite_N_Z_rows=n-finite,nonpositive_N_rows=nonpositive,h2_log_input_count=h2n,h2_log_count_matches_finite_rows=h2n==finite,historical_rg_template_count=rgcounts[0],rg_count_matches_total_rows=rgcounts[0]==n)
            rows.append(row);w.writerow(row);f.flush()
            print(f'COUNT_AUDIT {i}/100 {trait} template={n} finite={finite} elapsed={time.time()-started:.1f}',flush=True)
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'source_count':len(rows),'all100_expected_hashes_match':True,
        'all100_h2_nonmissing_counts_match':all(r['h2_log_count_matches_finite_rows'] for r in rows),
        'all1200_rg_template_counts_match':all(r['rg_count_matches_total_rows'] for r in rows),
        'code_sha256':sha(Path(__file__)),'table_sha256':sha(out),'elapsed_seconds':time.time()-started,
        'scientific_scope':'Processed input contents and stage-specific historical counts only; no native estimator or raw harmonization replay'}
    (PACKAGE/'logs/extension_processed_source_count_receipt_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
