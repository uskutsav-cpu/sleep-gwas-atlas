#!/usr/bin/env python3
"""Preserve originals; omit provenance-ineligible cells from scoring input.

No answer repair, hash replacement, model retry or frozen-instrument change.
Packet digest equality is a metadata gate, not proof of actual read access.
"""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--responses', type=Path, required=True)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--private-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    raw_bytes = args.responses.read_bytes()
    raw = [json.loads(x) for x in raw_bytes.decode().splitlines() if x.strip()]
    registration = json.loads(args.registration.read_text())
    schedule = {(x['model_alias'], x['condition'], x['task_id']): x
                for x in registration['schedule']}
    seen, eligible, excluded = set(), [], []
    for response in raw:
        ident = (response['model_alias'], response['condition'], response['task_id'])
        if ident in seen or ident not in schedule:
            raise ValueError('Duplicate/unregistered original cell')
        seen.add(ident)
        assigned = schedule[ident]
        if response['provider'] != 'OpenAI':
            raise ValueError('Incorrect provider metadata')
        if response.get('packet_sha256') != assigned['packet_sha256']:
            packet = args.private_root / assigned['packet_relative_path']
            actual = hashlib.sha256(packet.read_bytes()).hexdigest()
            if actual != assigned['packet_sha256']:
                raise ValueError('Assigned packet itself differs from registration')
            excluded.append(dict(
                model_alias=ident[0], condition=ident[1], task_id=ident[2],
                reason='UNEVALUABLE_PROVENANCE_RETURNED_PACKET_DIGEST_MISMATCH',
                returned_packet_sha256=response.get('packet_sha256'),
                registered_packet_sha256=assigned['packet_sha256'],
                assigned_file_sha256=actual,
                assigned_packet_relative_path=assigned['packet_relative_path'],
                original_response_relative_path=assigned['response_relative_path']))
        else:
            eligible.append(response)
    if seen != set(schedule):
        raise ValueError('Original definitive cell inventory incomplete')
    args.out.mkdir(parents=True, exist_ok=True)
    projected = args.out / 'eligible_responses.jsonl'
    projected.write_text(''.join(json.dumps(x, sort_keys=True) + '\n' for x in eligible))
    ledger = dict(
        original_definitive_cells=len(raw), eligible_cells=len(eligible),
        excluded_cells=len(excluded), exclusions=excluded,
        original_aggregate_sha256=hashlib.sha256(raw_bytes).hexdigest(),
        eligible_aggregate_sha256=hashlib.sha256(projected.read_bytes()).hexdigest(),
        registration_sha256=hashlib.sha256(args.registration.read_bytes()).hexdigest(),
        originals_modified=False, answers_repaired=False,
        metadata_gate_scope='Digest agreement only; no independent actual-access trace or technical blinding.',
        scoring_policy='Frozen scorer on unchanged eligible envelopes; excluded cells are missing/unevaluable provenance, not incorrect scientific answers.')
    (args.out / 'provenance_exclusion_ledger.json').write_text(json.dumps(ledger, indent=2) + '\n')
    print(json.dumps({k: ledger[k] for k in ['original_definitive_cells', 'eligible_cells', 'excluded_cells']}))


if __name__ == '__main__':
    main()
