#!/usr/bin/env python3
"""Recheck historical receipts and generate source metadata without native reruns.

Only package-owned tables/logs are written. Missing historical inputs are never
reconstructed under original filenames. Uses Python standard library only.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import re
import subprocess
from collections import Counter
import argparse
import time
from datetime import datetime, timezone

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parents[1]


def read(path):
    with (ROOT / path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write(name, rows):
    path = PACKAGE / 'tables' / name
    path.parent.mkdir(exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', restval='NA', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def effective_n(cases, controls, total):
    try:
        a, b = float(cases), float(controls)
        return format(4 / (1 / a + 1 / b), '.12g')
    except (ValueError, TypeError, ZeroDivisionError):
        return total if total not in ('PENDING', 'NOT_AVAILABLE') else 'NA'


def verify_remote_body(source_id, timeout):
    """Bounded current whole-body check; stream bytes to hashes, never to disk."""
    receipt = ROOT / 'discovery_extension/provenance/replication_streaming_receipts' / (source_id + '.json')
    original = json.loads(receipt.read_text())
    expected = original['source_verification']
    if expected['expected_size_bytes'] > 2 * 1024**3:
        raise ValueError('Single source exceeds the fixed 2 GiB verification cap')
    PACKAGE.joinpath('logs').mkdir(exist_ok=True)
    header_path = PACKAGE / 'logs' / ('provenance_remote_headers_' + source_id + '.txt')
    out = PACKAGE / 'logs' / ('provenance_remote_hash_' + source_id + '.json')
    if out.exists():
        prior = json.loads(out.read_text())
        suffix = re.sub(r'[^0-9]', '', prior['accessed_utc'])
        out.rename(out.with_name(out.stem + '_prior_' + suffix + '.json'))
        if header_path.exists():
            header_path.rename(header_path.with_name(header_path.stem + '_prior_' + suffix + '.txt'))
    # Account for both requested MVP objects and previous consumed attempts
    # before starting this campaign; no raw files are retained.
    target_ids = ['gwas_catalog_GCST90479148', 'gwas_catalog_GCST90479330']
    target_total = sum(json.loads((ROOT / 'discovery_extension/provenance/replication_streaming_receipts' / (sid + '.json')).read_text())['source_verification']['expected_size_bytes'] for sid in target_ids)
    prior_total = sum(json.loads(p.read_text())['bytes_streamed'] for p in (PACKAGE / 'logs').glob('provenance_remote_hash_*_prior_*.json'))
    if source_id in target_ids and target_total + prior_total > 2*1024**3:
        raise ValueError('Requested campaign plus prior consumed attempts exceeds 2 GiB')
    command = ['curl', '--silent', '--show-error', '--connect-timeout', '10', '--max-time', str(timeout),
               '--max-filesize', str(2 * 1024**3), '--dump-header', str(header_path), original['source_url']]
    started = time.monotonic(); next_progress = started + 10
    h = hashlib.sha256(); m = hashlib.md5(); count = 0
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        chunk = process.stdout.read(1048576)
        if not chunk:
            break
        count += len(chunk); h.update(chunk); m.update(chunk)
        if count > 2 * 1024**3:
            process.kill()
            break
        if time.monotonic() >= next_progress:
            print(json.dumps(dict(source_id=source_id, bytes_streamed=count, elapsed_seconds=round(time.monotonic()-started,2))), flush=True)
            next_progress = time.monotonic() + 10
    error = process.stderr.read().decode('utf-8', errors='replace'); code = process.wait()
    elapsed = time.monotonic() - started
    complete = code == 0 and count == expected['expected_size_bytes']
    match = complete and h.hexdigest() == expected['observed_sha256'] and m.hexdigest() == expected['observed_md5']
    result = dict(accessed_utc=datetime.now(timezone.utc).isoformat(), source_id=source_id,
                  historical_receipt=str(receipt.relative_to(ROOT)), source_url=original['source_url'],
                  bytes_streamed=count, elapsed_seconds=round(elapsed,3), curl_exit=code, error=error,
                  full_body_consumed=complete, sha256=h.hexdigest() if complete else None,
                  md5=m.hexdigest() if complete else None, partial_stream_sha256=None if complete else h.hexdigest(),
                  expected_sha256=expected['observed_sha256'], expected_md5=expected['observed_md5'],
                  current_raw_source_status='FULL_BODY_HASH_MATCH_HISTORICAL' if match else 'FULL_BODY_HASH_MISMATCH' if complete else 'INCOMPLETE_TIMEOUT_OR_ACCESS',
                  raw_body_retained=False, recovered_raw_source=False, native_gwas_rerun=False,
                  request_timeout_seconds=timeout, total_source_budget_bytes=2*1024**3,
                  two_mvp_expected_body_bytes=target_total, prior_attempt_bytes=prior_total,
                  response_header_file=str(header_path.relative_to(PACKAGE)))
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True), flush=True)


def main():
    checkpoint = json.loads((ROOT / 'discovery_extension/core_checkpoint.json').read_text())
    ledger = []

    def add(path, expected='NA', receipt='CURRENT_INVENTORY', evidence='PROCESSED_ARTIFACT', override=None):
        local = ROOT / path
        actual = digest(local) if local.is_file() else 'NA'
        status = ('CURRENT_FILE_HASHED' if expected == 'NA' else 'ORIGINAL_OUTPUT_HASH_VERIFIED' if actual == expected else 'HASH_MISMATCH') if local.is_file() else 'MISSING_ORIGINAL_SOURCE'
        ledger.append(dict(artifact_path=path, historical_expected_sha256=expected, current_sha256=actual,
                           status=override or status, evidence_level=evidence, receipt=receipt,
                           bytes=local.stat().st_size if local.is_file() else 'NA', native_rerun='NOT_RUN'))

    recovery = json.loads((PACKAGE / 'sources/provenance_recovery_scan.json').read_text())
    for path, expected in checkpoint['artifact_hashes_sha256'].items():
        recovered = [p for p in (PACKAGE / 'sources/recovered').glob('*') if expected in p.name]
        override = 'ORIGINAL_RECOVERED_IN_NEW_NAMESPACE' if recovered and not (ROOT / path).is_file() else None
        if path == 'environment/tool_versions.tsv' and (ROOT / path).is_file():
            old = subprocess.run(['git', '-C', str(ROOT), 'show', checkpoint['git_commit'] + ':' + path], text=True, capture_output=True, check=True)
            previous = list(csv.DictReader(io.StringIO(old.stdout), delimiter='\t'))
            current = read(path)
            ids = [r.get('component') for r in current]
            if current[:len(previous)] == previous and len(ids) == len(set(ids)):
                override = 'CHECKPOINT_APPEND_ONLY_ROWS_VERIFIED'
            else:
                override = 'CHECKPOINT_APPEND_ONLY_POLICY_FAIL'
        add(path, expected, 'discovery_extension/core_checkpoint.json', 'CORE_CHECKPOINT',
            override)
    final = json.loads((ROOT / 'discovery_extension/provenance/final_report.json').read_text())
    for group in ('inputs', 'outputs'):
        for entry in final[group].values():
            add(entry['path'], entry['sha256'] or 'NA', 'discovery_extension/provenance/final_report.json', 'FROZEN_EXTENSION',
                'ABSENT_BLOCKED_UPSTREAM' if entry['status'] == 'ABSENT_BLOCKED_UPSTREAM' else None)
    add('discovery_extension/results/ldsc/extension_pair_universe.tsv',
        json.loads((ROOT / 'discovery_extension/provenance/ldsc_collation.json').read_text())['pair_universe_sha256'],
        'discovery_extension/provenance/ldsc_collation.json')
    for path in sorted((ROOT / 'discovery_extension/provenance').rglob('*')):
        if path.is_file():
            add(str(path.relative_to(ROOT)), evidence='HISTORICAL_RECEIPT_CURRENT_BYTES_HASHED')
    for directory in ('streaming_receipts', 'replication_streaming_receipts'):
        for path in sorted((ROOT / 'discovery_extension/provenance' / directory).glob('*.json')):
            obj = json.loads(path.read_text()); receipt = str(path.relative_to(ROOT))
            for key in ('munged_output', 'harmonized_output', 'reference', 'hm3_alleles', 'harmonization_qc', 'munging_log', 'queue'):
                if key in obj and obj.get(key + '_sha256'):
                    add(obj[key], obj[key + '_sha256'], receipt, 'HISTORICAL_SOURCE_PROCESSING_RECEIPT')
            source = obj.get('source_verification', {})
            ledger.append(dict(artifact_path=obj.get('source_url', source.get('source', 'NA')),
                               historical_expected_sha256=source.get('observed_sha256', 'NA'), current_sha256='NA',
                               status='HISTORICAL_RECEIPT_ONLY_REMOTE_BODY_NOT_REDOWNLOADED', evidence_level='HISTORICAL_FULL_STREAM_RECEIPT',
                               receipt=receipt, bytes=source.get('observed_size_bytes', 'NA'), native_rerun='NOT_RUN'))
    for path in sorted((ROOT / 'discovery_extension/logs/replication').rglob('*.log')):
        add(str(path.relative_to(ROOT)), evidence='HISTORICAL_LDSC_LOG_CURRENT_BYTES_HASHED')
    for path in sorted((PACKAGE / 'sources/recovered').glob('*')):
        add(str(path.relative_to(ROOT)), path.name.split('_')[1], 'sources/provenance_recovery_scan.json', 'RECOVERED_ORIGINAL')
    remote_verification = {}
    for path in sorted((PACKAGE / 'logs').glob('provenance_remote_hash_*.json')):
        if '_prior_' in path.name:
            continue
        current = json.loads(path.read_text())
        remote_verification[current['source_id']] = current
        ledger.append(dict(artifact_path=current['source_url'], historical_expected_sha256=current['expected_sha256'],
                           current_sha256=current['sha256'] or 'NA', status='PASS_SOURCE_BODY_HASH' if current['current_raw_source_status']=='FULL_BODY_HASH_MATCH_HISTORICAL' else current['current_raw_source_status'],
                           evidence_level='CURRENT_REMOTE_BODY_VERIFICATION_NO_RAW_RETENTION',
                           receipt=str(path.relative_to(ROOT)), bytes=current['bytes_streamed'], native_rerun='BLOCKED_MISSING_GWAS_AND_LDSC_REFERENCES',
                           raw_body_retained='False'))
    # Preserve the append-only checkpoint policy even when this sparse worktree
    # has no environment file: compare original rows with authorized worktrees.
    snapshot = subprocess.run(['git', '-C', str(ROOT), 'show', checkpoint['git_commit'] + ':environment/tool_versions.tsv'], text=True, capture_output=True)
    env_checks = []
    if snapshot.returncode == 0:
        historical = list(csv.DictReader(io.StringIO(snapshot.stdout), delimiter='\t'))
        for directory in [str(ROOT)] + recovery['eligible_roots']:
            path = Path(directory) / 'environment/tool_versions.tsv'
            if not path.is_file():
                continue
            current = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'))
            ids = [row.get('component') for row in current]
            ok = current[:len(historical)] == historical and len(ids) == len(set(ids))
            by_id = {r.get('component'): r for r in current}
            env_checks.append(dict(path=str(path), current_sha256=digest(path), append_only_policy_pass=ok,
                                   original_row_values_unchanged=all(by_id.get(r.get('component'))==r for r in historical),
                                   historical_rows=len(historical), current_rows=len(current)))
    write('source_integrity_ledger.tsv', ledger)

    sleeps = [r for r in read('config/analysis_panel.tsv') if r['domain'] == 'sleep']
    modes = {r['trait_id']: r for r in read('config/sleep_measurement_modes.tsv')}
    schemas = {r['trait_id']: r for r in read('config/gwas_schemas.tsv')}
    sources = {r['source_id']: r for r in read('config/public_gwas_sources.tsv')}
    current_bibliography = {r['pmid']: r for r in json.loads((PACKAGE / 'sources/provenance_sleep_bibliography_corrections.json').read_text())['records']}
    for row in sleeps:
        row.update(modes[row['trait_id']])
        row['effective_n'] = effective_n(row['ncase'], row['ncontrol'], row['n_total'])
        row['cohort'] = 'FinnGen R9' if row['trait_id'] == 'sleep_apnea' else 'UK Biobank'
        row['effect_convention'] = schemas[row['trait_id']]['effect_convention']
        row['source_url'] = sources[row['source_id']]['download_url']
        row['historical_archive_sha256'] = sources[row['source_id']]['archive_sha256']
        row['current_dense_or_munged_source'] = 'UNAVAILABLE'
        row['metadata_evidence'] = 'LOCKED_PANEL_AND_SOURCE_SCHEMA; NOT_SOURCE_LEVEL_RERUN'
        row['verified_current_doi'] = current_bibliography.get(row['pmid'], {}).get('verified_doi', row['doi'])
        row['bibliography_verification_source'] = current_bibliography.get(row['pmid'], {}).get('primary_source', 'LOCKED_PANEL_METADATA')
    write('sleep_trait_metadata.tsv', sleeps)
    phenotypes = read('discovery_extension/config/candidate_traits.tsv')
    ph = {r['extension_trait_id']: r for r in phenotypes}
    queue = read('discovery_extension/results/replication/replication_source_queue.tsv')
    q = {r['pair_id']: r for r in queue}
    comparisons = []
    outcomes = read('discovery_extension/results/replication/replication_results.tsv')
    discovery = {(r['sleep_trait'], r['extension_trait_id']): r for r in read('discovery_extension/results/ldsc/extension_rg_matrix.tsv')}
    sm = {r['trait_id']: r for r in sleeps}
    mismatch_notes = {r['extension_trait_id']: r['phenotype_match_notes'] for p in
                      ('discovery_extension/config/replication_finngen_r13_sources.tsv', 'discovery_extension/config/replication_mvp_sources.tsv') for r in read(p)}
    for row in outcomes:
        entry = dict(row); source = q[row['pair_id']]; sleep = sm[row['sleep_trait']]; phenotype = ph[row['extension_trait_id']]
        d = discovery[(row['sleep_trait'], row['extension_trait_id'])]
        entry.update({k: source[k] for k in ('phenotype_match_status', 'replication_source_url', 'replication_checksum', 'replication_PMID',
                      'replication_DOI', 'sample_size', 'cases', 'controls', 'participant_overlap_evidence', 'effect_allele_status')})
        entry.update(sleep_source_id=sleep['source_id'], sleep_source_release=sleep['dataset_version'], sleep_source_cohort=sleep['cohort'],
                     sleep_source_pmid=sleep['pmid'], sleep_original_doi=sleep['doi'], sleep_verified_current_doi=sleep['verified_current_doi'],
                     sleep_phenotype_definition=sleep['phenotype_definition'], sleep_sample_size=sleep['n_total'], sleep_effective_n=sleep['effective_n'],
                     sleep_measurement_mode=sleep['collection_mode'], discovery_outcome_definition=phenotype['phenotype_definition'],
                     discovery_outcome_cases=phenotype['cases'], discovery_outcome_controls=phenotype['controls'], discovery_outcome_n=phenotype['sample_size'],
                     replication_effective_n=effective_n(source['cases'], source['controls'], source['sample_size']),
                     discovery_cross_trait_intercept=d['cross_trait_LDSC_intercept'], discovery_cross_trait_intercept_se=d['cross_trait_LDSC_intercept_se'],
                     discovery_valid_allele_snp_overlap=d['snp_overlap_valid_alleles'],
                     sleep_gwas_reused='YES_SAME_INPUT_PATH', sleep_input_path='data/munged/' + row['sleep_trait'] + '.sumstats.gz',
                     ld_reference='ref/eur_w_ld_chr/; historical European LD-score files; current bytes unavailable',
                     phenotype_comparison_notes=mismatch_notes.get(row['extension_trait_id'], 'NO_ELIGIBLE_SOURCE'),
                     independence_design='OUTCOME_COHORT_DISTINCT_SLEEP_GWAS_REUSED' if source['replication_source_id'] != 'NOT_AVAILABLE' else 'NOT_ESTIMABLE',
                     revised_positive_class='EXTERNAL_OUTCOME_SIDE_REPLICATION' if row['replication_class'] == 'REPLICATED' else 'NOT_A_POSITIVE_REPLICATION',
                     individual_level_overlap_verification='NOT_AVAILABLE_COHORT_DESIGN_EVIDENCE_ONLY',
                     estimate_covariance='UNKNOWN_SHARED_SLEEP_INPUT', native_source_rerun='BLOCKED_INPUTS_ABSENT')
        # A source can share FinnGen participants with R9 sleep apnea even if
        # it is independent of the UKB discovery outcome.
        if sleep['cohort'] == 'FinnGen R9' and source['replication_source_id'].startswith('finngen_'):
            entry['independence_design'] = 'PARTIAL_COHORT_OVERLAP_FINNGEN_R9_R13'
            if row['replication_class'] == 'REPLICATED':
                entry['revised_positive_class'] = 'PARTIAL_COHORT_OVERLAP'
        comparisons.append(entry)
    write('replication_source_comparison.tsv', comparisons)
    metadata = []
    for row in sleeps:
        metadata.append(dict(role='SLEEP', trait_id=row['trait_id'], source_id=row['source_id'], release=row['dataset_version'],
                             definition=row['phenotype_definition'], cohort=row['cohort'], ancestry=row['ancestry'], build=row['build'],
                             sample_size=row['n_total'], cases=row['ncase'], controls=row['ncontrol'], effective_n=row['effective_n'],
                             effect_convention=row['effect_convention'], source_url=row['source_url'], checksum=row['historical_archive_sha256'],
                             PMID=row['pmid'], original_DOI=row['doi'], verified_current_DOI=row['verified_current_doi'],
                             metadata_source='config/analysis_panel.tsv; config/public_gwas_sources.tsv'))
    for row in phenotypes:
        metadata.append(dict(role='DISCOVERY_OUTCOME', trait_id=row['extension_trait_id'], source_id=row['study_accession'], release=row['source'],
                             definition=row['phenotype_definition'], cohort='UK Biobank', ancestry=row['ancestry'], build=row['build'],
                             sample_size=row['sample_size'], cases=row['cases'], controls=row['controls'],
                             effective_n=effective_n(row['cases'], row['controls'], row['sample_size']), effect_convention='beta_EUR_per_alt',
                             source_url=row['source_url'], checksum=row['checksum'], metadata_source='discovery_extension/config/candidate_traits.tsv'))
        metadata[-1].update(PMID=row['PMID'], original_DOI=row['DOI'], verified_current_DOI=row['DOI'])
    unique = {r['replication_source_id']: r for r in queue if r['replication_source_id'] != 'NOT_AVAILABLE'}
    for source_id, row in unique.items():
        metadata.append(dict(role='VALIDATION_OUTCOME', trait_id=row['replication_study_accession'], source_id=source_id,
                             release='FinnGen R13' if source_id.startswith('finngen_') else 'MVP 2068-trait study',
                             definition=row['replication_phenotype_definition'], cohort='FinnGen' if source_id.startswith('finngen_') else 'MVP',
                             ancestry=row['ancestry'], build=row['build'], sample_size=row['sample_size'], cases=row['cases'], controls=row['controls'],
                             effective_n=effective_n(row['cases'], row['controls'], row['sample_size']), effect_convention='log_odds_per_effect_allele; FinnGen alt',
                             source_url=row['replication_source_url'], checksum=row['replication_checksum'],
                             metadata_source='discovery_extension/results/replication/replication_source_queue.tsv'))
        metadata[-1].update(PMID=row['replication_PMID'], original_DOI=row['replication_DOI'], verified_current_DOI=row['replication_DOI'])
        body_status = remote_verification.get(source_id, {}).get('current_raw_source_status', 'HISTORICAL_RECEIPT_ONLY')
        metadata[-1]['current_raw_source_integrity'] = 'PASS_SOURCE_BODY_HASH' if body_status=='FULL_BODY_HASH_MATCH_HISTORICAL' else body_status
    write('phenotype_and_source_metadata.tsv', metadata)
    original_missing_block = (ROOT / 'discovery_extension/analysis/frozen_extension_audit.md').read_text().split('Missing pinned core artifacts:')[1].split('Pinned core mismatches:')[0]
    formerly_missing = re.findall(r'`([^`]+)`', original_missing_block)
    newly_recovered = [path for path in formerly_missing if any(checkpoint['artifact_hashes_sha256'][path] in p.name for p in (PACKAGE / 'sources/recovered').glob('*'))]
    receipt = dict(schema_version='1.0.0', status='COMPLETED_WITH_EXTERNAL_SOURCE_BLOCKERS', native_gwas_reruns=0,
                   original_core_receipts=len(checkpoint['artifact_hashes_sha256']), source_integrity_rows=len(ledger),
                   source_integrity_status_counts=dict(Counter(r['status'] for r in ledger)),
                   historical_stream_receipts=sum(r['evidence_level']=='HISTORICAL_FULL_STREAM_RECEIPT' for r in ledger),
                   current_remote_body_statuses={sid:record['current_raw_source_status'] for sid,record in remote_verification.items()},
                   current_full_body_bytes=sum(record['bytes_streamed'] for record in remote_verification.values()),
                   current_raw_bodies_retained=sum(bool(record['raw_body_retained']) for record in remote_verification.values()),
                   environment_append_only_checks=env_checks, sleep_metadata_rows=len(sleeps), source_metadata_rows=len(metadata),
                   complete_replication_comparison_rows=len(comparisons),
                   positive_class_counts=dict(Counter(r['revised_positive_class'] for r in comparisons if r['replication_class']=='REPLICATED')),
                   measurement_positive_counts=dict(Counter(r['sleep_measurement_mode'] for r in comparisons if r['replication_class']=='REPLICATED')),
                   exact_original_missing_checkpoint_artifacts_recovered=len(newly_recovered), recovery_scan=recovery,
                   generated_table_hashes={name:digest(PACKAGE / 'tables' / name) for name in
                       ('source_integrity_ledger.tsv','phenotype_and_source_metadata.tsv','sleep_trait_metadata.tsv','replication_source_comparison.tsv')})
    (PACKAGE / 'logs').mkdir(exist_ok=True)
    (PACKAGE / 'logs/provenance_execution_receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('recovery_scan','environment_append_only_checks')}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-remote-body', choices=['gwas_catalog_GCST90479148', 'gwas_catalog_GCST90479330'], help='One of the two explicitly budgeted MVP sources; whole-body streaming without retaining raw bytes')
    parser.add_argument('--max-time', type=int, default=55, choices=range(1,3601), metavar='1..3600')
    args = parser.parse_args()
    if args.verify_remote_body:
        verify_remote_body(args.verify_remote_body, args.max_time)
    else:
        main()
