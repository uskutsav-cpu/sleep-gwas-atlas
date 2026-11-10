#!/usr/bin/env python3
"""Audit saved AlphaGenome per-variant files; no live model calls."""
import argparse, csv, gzip, hashlib, json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--queue', default='score_queue.csv')
p.add_argument('--store', default='alphagenome_scores')
p.add_argument('--output', default='score_coverage_audit.csv')
a = p.parse_args()
queue = list(csv.DictReader(open(a.queue, newline='')))
if not queue or len({x['variant_key'] for x in queue}) != len(queue):
    raise SystemExit('Queue empty or variant keys are duplicated')
store = Path(a.store)
config_path = store / 'run_config.json'
config = json.loads(config_path.read_text()) if config_path.exists() else None
fields = ('client_version','model_version','organism','genome_build','sequence_length','scorers')
expected = set(config.get('scorers', [])) if config else set()
results = []
for item in queue:
    key = item['variant_key']
    digest = hashlib.sha256(key.encode()).hexdigest()[:20]
    receipt_path, score_path = store / (digest+'.json'), store / (digest+'.csv.gz')
    result = {'variant_key': key, 'status': 'missing', 'rows': 0, 'returned_scorers': '', 'returned_output_types': '', 'missing_output_types': '', 'detail': ''}
    if receipt_path.exists():
        try:
            receipt = json.loads(receipt_path.read_text())
            if receipt.get('variant_key') != key:
                raise ValueError('receipt variant key mismatch')
            if receipt.get('status') != 'complete':
                result['status'] = 'failed_or_incomplete'
            elif not config:
                raise ValueError('run configuration missing')
            elif any(receipt.get('config',{}).get(k) != config.get(k) for k in fields):
                raise ValueError('receipt configuration differs from run configuration')
            elif not score_path.is_file():
                raise ValueError('score file missing')
            elif hashlib.sha256(score_path.read_bytes()).hexdigest() != receipt.get('sha256'):
                raise ValueError('score file hash mismatch')
            else:
                with gzip.open(score_path,'rt',newline='') as f:
                    reader = csv.DictReader(f)
                    required = {'variant_id','output_type','variant_scorer','raw_score'}
                    if not required.issubset(reader.fieldnames or []):
                        raise ValueError('required score columns missing')
                    scorers, output_types, ids, count = set(), set(), set(), 0
                    for row in reader:
                        count += 1
                        if row['variant_scorer']: scorers.add(row['variant_scorer'])
                        if row['output_type']: output_types.add(row['output_type'])
                        ids.add(row['variant_id'])
                if not count or count != receipt.get('rows'):
                    raise ValueError('empty file or receipt row count mismatch')
                if ids != {key}:
                    raise ValueError('score variant identifiers do not match queue')
                result.update(status='verified_file',rows=count,returned_scorers='|'.join(sorted(scorers)),returned_output_types='|'.join(sorted(output_types)),missing_output_types='|'.join(sorted(expected-output_types)))
                if not expected:
                    result['detail']='No scorer names in configuration; coverage unavailable'
                elif expected-output_types:
                    result['detail']='At least one requested output type was not returned'
        except (OSError, ValueError, TypeError, KeyError, EOFError, json.JSONDecodeError) as e:
            result.update(status='invalid', detail=str(e))
    results.append(result)
with open(a.output,'w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
for state in ('verified_file','missing','failed_or_incomplete','invalid'):
    print(state, sum(r['status']==state for r in results))
print('Output-type coverage does not verify exact scorer settings or track completeness; tissue, gene, event and biological validity are not audited.')
