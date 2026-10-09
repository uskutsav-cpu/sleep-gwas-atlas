#!/usr/bin/env python3
"""Read-only SSD discovery and byte-identical checkpoint recovery.

Only output_dir is writable. No input is renamed, repaired, or deleted.
Raw GWAS files are never copied into the repository.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from datetime import datetime, timezone

SSD = Path('/Volumes/Extreme SSD')
ARCHIVE = SSD / 'Utsav-Research-Archive/Sleep-GWAS'
FAILED = ARCHIVE / 'FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
DEPS = ARCHIVE / 'Sep1-local-dependencies'
RAW = SSD / 'Codex Archive/2026-08-26-sleep-gwas-atlas'
PRIOR = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/goa/work/sleep-gwas-atlas')
OLD = Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas')

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def command(args):
    p = subprocess.run(args, capture_output=True, text=True)
    return {'command': args, 'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    repo, out = args.repo.resolve(), args.output.resolve()
    for sub in ['logs', 'manifests', 'hashes', 'tables', 'figures', 'tests', 'sources/recovered']:
        (out / sub).mkdir(parents=True, exist_ok=True)
    checkpoint = json.loads((repo / 'discovery_extension/core_checkpoint.json').read_text())
    audit, recovery = [], []
    for relative, expected in checkpoint['artifact_hashes_sha256'].items():
        candidates = [FAILED / relative, PRIOR / relative, OLD / relative, repo / relative]
        prior_recovered = PRIOR / 'discovery_extension/sleep_submission_evidence_v1/sources/recovered'
        candidates += sorted(prior_recovered.glob(f'core_{expected}_*'))
        match = None
        for candidate in candidates:
            if not candidate.is_file():
                continue
            actual = sha256(candidate)
            row = {'artifact': relative, 'path': str(candidate), 'bytes': candidate.stat().st_size,
                   'expected_sha256': expected, 'actual_sha256': actual,
                   'status': 'EXACT' if actual == expected else 'DIFFERENT_BYTES'}
            audit.append(row)
            if actual == expected and match is None:
                match = candidate
        record = {'artifact': relative, 'expected_sha256': expected, 'status': 'UNRECOVERED',
                  'source_path': '', 'actual_sha256': '', 'recovered_path': '', 'bytes': '',
                  'recomputed': 'false'}
        if match is not None:
            destination = out / 'sources/recovered' / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and sha256(destination) != expected:
                raise RuntimeError(f'Refusing to overwrite differing recovered bytes: {destination}')
            if not destination.exists():
                shutil.copyfile(match, destination)
            assert sha256(destination) == expected
            record.update(status='EXACT_SHA256_RECOVERED', source_path=str(match),
                          actual_sha256=expected, recovered_path=str(destination.relative_to(out)),
                          bytes=destination.stat().st_size)
        recovery.append(record)
    write_tsv(out / 'ORIGINAL_CORE_RECOVERY.tsv', recovery, list(recovery[0]))
    write_tsv(out / 'tables/core_candidate_hash_audit.tsv', audit, list(audit[0]))

    panel = read_tsv(repo / 'config/analysis_panel.tsv')
    schema = {r['trait_id']: r for r in read_tsv(repo / 'config/gwas_schemas.tsv')}
    sources = read_tsv(repo / 'config/public_gwas_sources.tsv')
    raw_traits = {r['raw_file']: r for r in panel}
    archive_sources = {r['archive_name']: r for r in sources}
    inventory, errors = [], []
    roots = [RAW, DEPS / 'data', FAILED / 'results', FAILED / 'config', FAILED / 'ref',
             FAILED / 'reference', FAILED / 'environment', FAILED / 'discovery_extension',
             ARCHIVE / 'Brain6-Replication-R13/hm3', ARCHIVE / 'Brain6-Replication-R13/munged']
    for root in roots:
        if not root.exists():
            errors.append({'path': str(root), 'status': 'ABSENT_NAMESPACE'})
            continue
        def onerror(error):
            errors.append({'path': error.filename, 'status': str(error)})
        for directory, dirs, files in os.walk(root, followlinks=False, onerror=onerror):
            dirs[:] = sorted(d for d in dirs if d not in ['track_b', '_smoketest', '_test_raw', '_test_harmonized', '.git'] and not d.startswith('._'))
            for name in sorted(files):
                if name.startswith('._') or name in ['.DS_Store', '.gitkeep']:
                    continue
                path = Path(directory) / name
                if path.is_symlink():
                    errors.append({'path': str(path), 'status': 'SYMLINK_NOT_FOLLOWED'})
                    continue
                stat = path.stat()
                trait = raw_traits.get(name, {})
                if name.endswith('.sumstats.gz'):
                    trait = next((r for r in panel if name == r['trait_id'] + '.sumstats.gz'), {})
                if name.endswith('.harmonized.tsv.gz'):
                    trait = next((r for r in panel if name.startswith(r['trait_id'] + '.')), {})
                source = archive_sources.get(name, {})
                expected = source.get('archive_sha256', '')
                kind = 'auxiliary_or_historical_output'
                if '/raw/' in str(path): kind = 'raw_gwas'
                if '/harmonized/' in str(path): kind = 'harmonized_gwas'
                if '/munged/' in str(path) and name.endswith('.sumstats.gz'): kind = 'hapmap3_munged_gwas'
                if '/munged_chromosome_split/' in str(path): kind = 'chromosome_munged_gwas'
                actual = ''
                if stat.st_size <= 1024 * 1024 and kind == 'auxiliary_or_historical_output':
                    actual = sha256(path)
                status = 'CURRENT_HASH_NO_EXPECTED_HASH' if actual else 'INVENTORIED_NOT_YET_HASHED'
                inventory.append({'path': str(path), 'bytes': stat.st_size,
                    'expected_sha256': expected, 'actual_sha256': actual,
                    'source_release': trait.get('dataset_version', source.get('source_id', 'UNKNOWN')),
                    'genome_build': trait.get('build', source.get('build_status', 'UNKNOWN')),
                    'ancestry': trait.get('ancestry', source.get('ancestry_reported', 'UNKNOWN')),
                    'phenotype_identity': trait.get('trait_id', source.get('trait_ids', 'UNKNOWN')),
                    'effect_coding': schema.get(trait.get('trait_id', ''), {}).get('effect_convention', 'UNKNOWN'),
                    'permitted_use': 'RESEARCH_ONLY; SOURCE_LICENSE_REVIEW_REQUIRED; NO_RAW_REDISTRIBUTION',
                    'required_downstream_tasks': 'VERIFY_SOURCE_AND_PREPROCESSING_PROVENANCE',
                    'verification_status': status, 'asset_kind': kind})
    fields = ['path','bytes','expected_sha256','actual_sha256','source_release','genome_build','ancestry',
              'phenotype_identity','effect_coding','permitted_use','required_downstream_tasks','verification_status','asset_kind']
    write_tsv(out / 'SSD_INPUT_MANIFEST.tsv', inventory, fields)
    (out / 'logs/inventory_errors.json').write_text(json.dumps(errors, indent=2) + '\n')

    commands = [['diskutil','info',str(SSD)], ['df','-k',str(SSD),str(repo)],
                ['sysctl','hw.memsize','hw.ncpu','vm.swapusage'], ['vm_stat'],
                ['git','-C',str(repo),'branch','-avv'], ['git','-C',str(OLD),'worktree','list','--porcelain'],
                ['git','-C',str(OLD),'status','--short'], ['git','--version'], ['gh','--version']]
    records = [command(c) for c in commands]
    processes = command(['ps','-axo','pid,ppid,%cpu,%mem,etime,command'])
    processes['stdout'] = '\n'.join(s for s in processes['stdout'].splitlines()
        if any(t in s.lower() for t in ['ldsc', 'munge_sumstats', 'sleep-gwas', 'genomicsem', 'rscrip', 'snakemake']))
    records.append(processes)
    records.append(command(['lsof','-nP',str(SSD)]))
    (out / 'logs/live_recovery_preflight.json').write_text(json.dumps(records, indent=2) + '\n')

    # Compare the immutable environment row order explicitly; do not repair historical files.
    env = out / 'sources/recovered/environment/tool_versions.tsv'
    if env.exists():
        original = read_tsv(env)
        comparisons = []
        for base in [FAILED, PRIOR, OLD, repo]:
            current = base / 'environment/tool_versions.tsv'
            if not current.exists(): continue
            rows = read_tsv(current)
            comparisons.append({'path': str(current), 'sha256': sha256(current),
                'exact_bytes': sha256(current) == checkpoint['artifact_hashes_sha256']['environment/tool_versions.tsv'],
                'original_rows_present_with_unchanged_values': all(r in rows for r in original),
                'immutable_prefix_matches': rows[:len(original)] == original,
                'original_rows': len(original), 'current_rows': len(rows)})
        (out / 'manifests/environment_order_audit_v1.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    summary = {'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'ssd_mount': str(SSD), 'input_namespace_roots': list(map(str, roots)),
        'inventory_records': len(inventory), 'inventory_bytes': sum(r['bytes'] for r in inventory),
        'checkpoint_artifacts': len(recovery), 'exact_recovered': sum(r['status']=='EXACT_SHA256_RECOVERED' for r in recovery),
        'unrecovered': [r['artifact'] for r in recovery if r['status']=='UNRECOVERED'],
        'input_files_modified': False, 'hash_verification_not_native_scientific_reproduction': True}
    (out / 'logs/recovery_receipt_v1.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
