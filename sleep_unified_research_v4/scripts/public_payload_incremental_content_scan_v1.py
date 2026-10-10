"""Bounded credential/body-name check of the current Git-selected research payload.

Immutable scan receipts are curation evidence, not source-rights clearance.
A prior receipt permits the exact same content identities to reuse their prior
pattern result; new or changed content is scanned. No file is staged or sent.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import resource
import subprocess
import time
from datetime import datetime, timezone

PATTERNS = {
    'private_key_header': r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
    'github_token_literal': r'\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{30,}\b',
    'aws_access_key_literal': r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
    'openai_key_literal': r'\bsk-(?:proj-)?[A-Za-z0-9_-]{35,}\b',
}
TEXT_SUFFIXES = {'.md', '.json', '.log', '.tsv', '.py', '.diff', '.txt', '.sha256', '.jsonl', ''}
FORBIDDEN_NAME_PARTS = ('.sumstats', '.vcf', '.bgen', '.bed', '.bim', '.fam', '.bcor', '.panel', '.parquet', '.pickle', '.npy', '.npz')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--prior', type=Path)
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    root = package.parent
    output = args.output.resolve()
    assert output.parent == package / 'logs'
    assert not output.exists() and not output.is_symlink()
    prior = {}
    if args.prior:
        prior_body = args.prior.read_bytes()
        prior_receipt = json.loads(prior_body)
        assert prior_receipt['credential_literal_findings'] == []
        assert prior_receipt['credential_literal_pattern_definitions'] == PATTERNS
        prior = prior_receipt['file_identity']
    started = time.monotonic()
    names = subprocess.run(
        ['git', 'ls-files', '--others', '--exclude-standard', '--', package.name + '/'],
        cwd=root, text=True, capture_output=True, check=True,
    ).stdout.splitlines()
    # Relevant modifications already tracked before this research package.
    names += subprocess.run(
        ['git', 'diff', '--name-only', '--', package.name + '/', '.gitignore', 'README.md'],
        cwd=root, text=True, capture_output=True, check=True,
    ).stdout.splitlines()
    files = {}
    findings = []
    forbidden = []
    extra = []
    reused = 0
    rescanned = 0
    for name in sorted(set(names)):
        q = root / name
        assert q.is_file() and not q.is_symlink()
        assert not any(parent.is_symlink() for parent in q.parents)
        assert q.stat().st_size <= 20 << 20, name
        body = q.read_bytes()
        digest = hashlib.sha256(body).hexdigest()
        row = {'bytes': len(body), 'sha256': digest,
               'kind': 'text' if q.suffix in TEXT_SUFFIXES else 'scientific_figure_binary'}
        files[name] = row
        if any(token in q.name.lower() for token in FORBIDDEN_NAME_PARTS):
            forbidden.append(name)
        if q.suffix in TEXT_SUFFIXES:
            if prior.get(name) == row:
                reused += 1
            else:
                s = body.decode('utf-8', errors='strict')
                for kind, pattern in PATTERNS.items():
                    for match in re.finditer(pattern, s):
                        findings.append({'path': name, 'kind': kind,
                                         'line': s[:match.start()].count('\n') + 1})
                rescanned += 1
        if q.suffix in ('', '.txt', '.jsonl'):
            extra.append({'path': name, 'bytes': len(body)})
    receipt = {
        'schema': 'preliminary_incremental_public_payload_content_scan_v1',
        'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS_BOUNDED_CONTENT_PATTERNS' if not findings and not forbidden else 'STOPPED_NO_PUBLICATION',
        'candidate_file_count': len(files),
        'candidate_regular_bytes': sum(v['bytes'] for v in files.values()),
        'file_identity': files,
        'credential_literal_pattern_definitions': PATTERNS,
        'credential_literal_findings': findings,
        'forbidden_raw_genotype_or_runtime_body_filename_findings': forbidden,
        'symlink_candidates': 0,
        'text_files_reusing_exact_prior_pattern_result': reused,
        'new_or_changed_text_files_scanned': rescanned,
        'files_with_unusual_text_suffix': extra,
        'prior_receipt_sha256': hashlib.sha256(prior_body).hexdigest() if args.prior else None,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'current_scope': 'Git-selected untracked research_v4 files plus relevant current tracked modifications only. Future artifacts, unchanged tracked repository content, redistribution rights and final publication choice are separate gates.',
        'checks_do_not_prove_absence_of_all_secrets_or_participant_information': True,
        'source_permissions_or_human_review_cleared': False,
        'Git_staged_committed_or_pushed': False,
        'elapsed_seconds': time.monotonic() - started,
        'peak_parent_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    with output.open('x') as stream:
        stream.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'candidate_file_count': len(files),
                      'receipt_path': str(output),
                      'receipt_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                      'credential_findings': findings, 'forbidden_filename_findings': forbidden,
                      'new_or_changed_text_files_scanned': rescanned, 'reused_text_checks': reused,
                      'unusual_text_suffix_files': extra}))
    if findings or forbidden:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
