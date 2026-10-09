#!/usr/bin/env python3
"""Bounded read-only archive scan; retain only exact missing source hashes."""
from collections import defaultdict
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import time

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
V1=ROOT/'sleep_unified_research_v1'
ARCHIVE=Path('/Volumes/Extreme SSD/Codex-Archive/2026-08-26.tar.gz')
DEST=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/missing_core_sources_v2')
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()

def main():
    for sub in ['manifests','tables','logs']: (PACKAGE/sub).mkdir(parents=True,exist_ok=True)
    receipt=PACKAGE/'logs/missing_core_archive_search_receipt_v2.json'
    if receipt.exists():raise SystemExit('PRIOR_SEARCH_PRESERVED')
    missing=[r for r in read(V1/'tables/native_input_hash_checks.tsv') if r['status']=='MISSING']
    sources=read(ROOT/'config/public_gwas_sources.tsv')
    source_by_sha={r['archive_sha256']:r for r in sources if r['archive_sha256']}
    targets=defaultdict(list);sizes={}
    for r in missing:
        targets[r['expected_sha256']].append(r)
        s=source_by_sha.get(r['expected_sha256'])
        if s:sizes[r['expected_sha256']]=int(s['archive_bytes'])
    by_name={Path(r['path']).name:h for h,rs in targets.items() for r in rs}
    DEST.mkdir(parents=True,exist_ok=True)
    budget=sum(sizes.values())+1024**3
    if shutil.disk_usage(DEST).free<budget+5*1024**3:raise SystemExit('SSD_RESOURCE_GATE_FAILED')
    plan={'frozen_utc':datetime.now(timezone.utc).isoformat(),'archive':str(ARCHIVE),'archive_bytes':ARCHIVE.stat().st_size,
          'missing_ledger_entries':len(missing),'unique_expected_hashes':len(targets),
          'expected_unique_source_bytes':sum(sizes.values()),'retention_budget_bytes':budget,
          'expected_sizes':sizes,'target_names':by_name,'targets_by_sha256':dict(targets),
          'reader':'one sequential tar/gzip reader; <=4MB buffers; member metadata discarded after use',
          'network_bytes':0,'destination':str(DEST),'original_SSD_and_v1_outputs_modified':False}
    plan_path=PACKAGE/'manifests/missing_core_archive_search_plan_v2.json'
    plan_path.write_text(json.dumps(plan,indent=2)+'\n')
    print(json.dumps({'plan':str(plan_path),'targets':len(targets),'retention_budget_bytes':budget}),flush=True)
    rows=[];found={};started=time.time();member_count=0
    table=PACKAGE/'tables/missing_core_archive_candidates_v2.tsv'
    if table.exists():raise SystemExit('PRIOR_PARTIAL_SEARCH_PRESERVED')
    fields=['archive_member','type','bytes','target_sha256','actual_sha256','status','retained_path','alias_entries_satisfied']
    with table.open('x',newline='') as f,tarfile.open(ARCHIVE,'r|gz') as archive:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for member in archive:
            member_count+=1
            # TarFile's streaming reader otherwise accumulates all TarInfo records.
            archive.members.clear()
            name=Path(member.name).name
            if name.startswith('._') or name not in by_name:continue
            expected=by_name[name]
            row=dict(archive_member=member.name,type=member.type.decode('ascii',errors='replace'),bytes=member.size,target_sha256=expected,actual_sha256='',status='',retained_path='',alias_entries_satisfied='')
            if not member.isfile():row['status']='NONREGULAR_MEMBER_NOT_EXTRACTED'
            elif expected in sizes and member.size!=sizes[expected]:row['status']='DIFFERENT_BYTE_COUNT_CANNOT_MATCH'
            elif expected in found:
                row.update(status='TARGET_ALREADY_RECOVERED; DUPLICATE_MEMBER_NOT_REHASHED',retained_path=found[expected]['retained_path'])
            else:
                temporary=DEST/(expected[:16]+'.extracting')
                if temporary.exists():raise SystemExit('PRIOR_PARTIAL_SOURCE_PRESERVED')
                h=hashlib.sha256();n=0
                with archive.extractfile(member) as source,temporary.open('xb') as destination:
                    for b in iter(lambda:source.read(4*1024**2),b''):
                        h.update(b);destination.write(b);n+=len(b)
                row['actual_sha256']=h.hexdigest()
                if n!=member.size:raise SystemExit('TRUNCATED_MEMBER_PRESERVED')
                if row['actual_sha256']==expected:
                    final=DEST/(expected[:16]+'__'+name)
                    if final.exists():
                        if sha(final)!=expected:raise SystemExit('EXISTING_RECOVERY_DIFFERS; NO_OVERWRITE')
                        temporary.unlink()
                    else:temporary.rename(final)
                    row.update(status='EXACT_EXPECTED_SHA256_RECOVERED',retained_path=str(final),alias_entries_satisfied=len(targets[expected]));found[expected]=row
                else:
                    # Only this task-created failed extraction is removed; original inputs stay untouched.
                    temporary.unlink();row['status']='DIFFERENT_SHA256_NOT_RETAINED'
            rows.append(row);w.writerow(row);f.flush()
            print(f"ARCHIVE_CANDIDATE {member_count} {name} {row['status']} elapsed={time.time()-started:.1f}",flush=True)
    completed={'completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.time()-started,
        'archive':str(ARCHIVE),'archive_members_scanned':member_count,'candidate_members':len(rows),
        'unique_sources_recovered':len(found),'ledger_entries_satisfied':sum(len(targets[h]) for h in found),
        'remaining_expected_hashes':[h for h in targets if h not in found],
        'recovered_sources':found,'table_sha256':sha(table),'plan_sha256':sha(plan_path),'script_sha256':sha(Path(__file__)),
        'original_SSD_and_v1_outputs_modified':False,'native_preprocessing_or_estimation_executed':False}
    receipt.write_text(json.dumps(completed,indent=2)+'\n')
    print(json.dumps(completed,indent=2),flush=True)
if __name__=='__main__':main()
