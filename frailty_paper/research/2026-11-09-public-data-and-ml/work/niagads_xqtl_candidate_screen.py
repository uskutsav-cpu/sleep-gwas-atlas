#!/usr/bin/env python3
"""Retrieve publicly exposed significant association records by candidate rsID."""
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent.parent
out = root / 'outputs'
files = [out / 'unvalidated_15_variant_score_prototype.csv', out / 'unvalidated_22_shared_direction_score.csv']
ids = set()
for path in files:
    with path.open(newline='') as f:
        ids.update(row['variant'] for row in csv.DictReader(f))
assert len(ids) == 33

results = {}
for rsid in ['rs6733839'] + sorted(ids):
    url = 'https://xqtl.niagads.org/download/rsid/' + rsid
    try:
        result = subprocess.run(['curl', '--fail', '--silent', '--show-error', '--max-time', '30', url],
                                capture_output=True, check=True, timeout=35)
        body = result.stdout.decode('utf-8')
        if 'No records found' in body:
            results[rsid] = {'status': 'no_hmt_records_returned', 'rows': []}
        elif body.startswith('rsID\tvariantId\t'):
            rows = list(csv.DictReader(io.StringIO(body), delimiter='\t'))
            results[rsid] = {'status': 'records_returned', 'rows': rows}
        else:
            results[rsid] = {'status': 'unexpected_response', 'prefix': body[:200]}
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, UnicodeError) as exc:
        results[rsid] = {'status': 'request_failed', 'error': str(exc)}
receipt = {'release': 'NIAGADS ADSP FunGen xQTL Atlas v1.0, portal HMT-filtered',
           'query_time_utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
           'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
           'positive_control': 'rs6733839', 'results': results,
           'caveat': 'Returned rows are selected associations; no result is not a null test. rsID matching alone does not establish the exact project allelic contrast.'}
(out / 'stage6_niagads_xqtl_candidate_screen.json').write_text(json.dumps(receipt, indent=2) + '\n')
for name, result in results.items():
    print(name, result['status'], len(result.get('rows', [])))
