#!/usr/bin/env python3
"""Recover only exact receipt-pinned native inputs into a new SSD directory.

Uses a streaming tar reader; no archive path is trusted as an output path.
Source archives and existing SSD files are read-only. Input summaries stay on
the SSD and never enter the Git commit.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile
from datetime import datetime, timezone

PREFIX = '2026-08-26/go/work/sleep-gwas-atlas/'
DEFAULT_ARCHIVE = Path('/Volumes/Extreme SSD/Codex-Archive/2026-08-26.tar.gz')
DEFAULT_DEST = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''): h.update(block)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument('--archive', type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument('--destination', type=Path, default=DEFAULT_DEST)
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args(); repo, dest = args.repo.resolve(), args.destination.resolve()
    out = repo/'sleep_unified_research_v1'
    expected, evidence = {}, {}
    for namespace in ['streaming_receipts','replication_streaming_receipts']:
        for receipt in sorted((repo/'discovery_extension/provenance'/namespace).glob('*.json')):
            r = json.loads(receipt.read_text())
            for path_key, hash_key in [('munged_output','munged_output_sha256'),('munging_log','munging_log_sha256'),('harmonization_qc','harmonization_qc_sha256')]:
                rel, sha = r.get(path_key),r.get(hash_key)
                if rel and sha:
                    relative = Path(rel)
                    if relative.is_absolute() or '..' in relative.parts or relative.parts[0]!='discovery_extension':
                        raise RuntimeError('Unsafe historical receipt path')
                    if relative.as_posix() in expected and expected[relative.as_posix()] != sha:
                        raise RuntimeError('Conflicting expected hashes')
                    expected[relative.as_posix()] = sha
                    evidence[relative.as_posix()] = str(receipt.relative_to(repo))
    reference = json.loads((repo/'discovery_extension/provenance/panukbb/hm3_variant_reference_build.json').read_text())
    expected[reference['output']] = reference['output_sha256']
    evidence[reference['output']] = 'discovery_extension/provenance/panukbb/hm3_variant_reference_build.json'
    # Native original logs are recovered as independent audit assets. Their
    # hashes are compared below only when existing historical receipts pin them.
    extra_patterns = ['discovery_extension/logs/h2/','discovery_extension/logs/rg/',
                      'discovery_extension/logs/replication/h2/', 'discovery_extension/logs/replication/rg/']
    plan = {'archive':str(args.archive),'archive_bytes':args.archive.stat().st_size,
            'destination':str(dest),'expected_receipt_pinned_artifacts':len(expected),
            'disk_retention_plan':'Only ~1.2 GB of munged inputs/reference/logs/QC; no full dense raw extraction.',
            'network_transfer_bytes':0,'reader_mode':'one sequential gzip/tar pass; one writer; <=4MB buffers',
            'original_ssd_inputs_modified':False,'large_data_in_git':False}
    (out/'manifests/archive_extraction_plan_v1.json').write_text(json.dumps(plan,indent=2)+'\n')
    if args.dry_run:
        print(json.dumps(plan,indent=2)); return
    if not str(dest).startswith('/Volumes/Extreme SSD/sleep-unified-research-v1/'):
        raise SystemExit('Destination must be in this task\'s new SSD namespace')
    dest.mkdir(parents=True,exist_ok=True)
    rows = []; recovered = set()
    with tarfile.open(args.archive, 'r|gz') as tar:
        for member in tar:
            if not member.isfile() or not member.name.startswith(PREFIX): continue
            rel = member.name[len(PREFIX):]
            relative = Path(rel)
            if relative.is_absolute() or '..' in relative.parts or relative.name.startswith('._'): continue
            if rel not in expected and not (any(rel.startswith(x) for x in extra_patterns) and rel.endswith('.log')): continue
            output = dest/relative
            output.parent.mkdir(parents=True,exist_ok=True)
            # A previous extraction may resume only with exact expected bytes.
            if output.exists():
                actual = digest(output)
                if expected.get(rel) and actual != expected[rel]:
                    raise RuntimeError('Refusing to replace differing recovered bytes')
            else:
                temporary = output.with_name(output.name+'.extracting')
                if temporary.exists():
                    raise RuntimeError('Prior partial extraction requires explicit inspection')
                h=hashlib.sha256()
                with tar.extractfile(member) as src, temporary.open('xb') as dst:
                    for block in iter(lambda:src.read(4*1024*1024),b''):
                        h.update(block); dst.write(block)
                actual=h.hexdigest()
                if expected.get(rel) and actual != expected[rel]:
                    raise RuntimeError(f'ARCHIVED_INPUT_HASH_MISMATCH {rel}: partial retained as {temporary}')
                temporary.rename(output)
            row={'archive_path':str(args.archive),'archive_member':member.name,'path':str(output),
                 'bytes':member.size,'expected_sha256':expected.get(rel,''),'actual_sha256':actual,
                 'expected_hash_evidence':evidence.get(rel,''),
                 'status':'EXACT_RECEIPT_HASH_RECOVERED' if rel in expected else 'HISTORICAL_LOG_RECOVERED_CURRENT_HASH_ONLY'}
            rows.append(row); recovered.add(rel)
            with (out/'tables/archived_extension_recovery.tsv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t'); w.writeheader(); w.writerows(rows)
            print(f"RECOVERED {len(rows)} {rel} {row['status']}",flush=True)
    missing=sorted(set(expected)-recovered)
    result={**plan,'completed_utc':datetime.now(timezone.utc).isoformat(),'recovered_artifacts':len(rows),
            'recovered_receipt_pinned_artifacts':len(set(expected)&recovered),
            'missing_receipt_pinned_artifacts':missing,'recovered_bytes':sum(r['bytes'] for r in rows),
            'extension_munged_sources':sum('/data/munged/' in r['path'] and r['path'].endswith('.sumstats.gz') for r in rows),
            'replication_munged_sources':sum('/data/replication/munged/' in r['path'] and r['path'].endswith('.sumstats.gz') for r in rows)}
    (out/'logs/archived_extension_recovery_receipt_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
