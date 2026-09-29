#!/usr/bin/env python3
"""Read-only, non-admission schema crosswalk for archived protected Track B."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RECOVERY = ROOT / 'brain6/results/track_b_recovery_v1/track_b_recovery_audit.json'
TEMPLATE = ROOT / 'brain6/results/placo/insomnia__mdd/variants.tsv.gz'
LOCK = ROOT / 'brain6/config/placo_family_v3/family_lock.json'
SOURCE_HEADER = ['analysis_id','pair_id','SNP','CHR','BP','A1','A2','Z1','Z2','P1','P2',
                 'T_PLACO_PLUS','P_PLACO_PLUS','PLACO_BH_Q','within_pair_family_n',
                 'analysis_status','numerical_error']
TARGET_HEADER = ['SNP','CHR','BP','Z1','Z2','P_PLACO','Q_WITHIN_PAIR',
                 'Q_BONFERRONI_ACROSS_PAIRS','status','headline']


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def checked_float(raw: str, field: str, row: int, probability: bool = False) -> float:
    value = float(raw)
    if not math.isfinite(value) or (probability and not 0 <= value <= 1):
        raise ValueError(f'invalid {field} at row {row}')
    return value


def stable_write(path: Path, data: str) -> None:
    if path.exists() and path.read_text() != data:
        raise RuntimeError(f'existing output differs: {path}')
    if not path.exists():
        path.write_text(data)


def main() -> None:
    recovery = json.loads(RECOVERY.read_text())
    family = json.loads(LOCK.read_text())
    ledger = Path(recovery['file_hashes']['archived_B_ledger']['path'])
    expected_sha = recovery['file_hashes']['archived_B_ledger']['sha256']
    if sha256(ledger) != expected_sha:
        raise ValueError('archived B ledger checksum changed')
    with gzip.open(TEMPLATE, 'rt', newline='') as f:
        actual_target_header = next(csv.reader(f, delimiter='\t'))
    if actual_target_header != TARGET_HEADER or family['protected_legacy_pair'] != 'insomnia__adhd':
        raise ValueError('v3 target schema or family lock changed')
    if family['n_selected_tracks'] != 5 or family['headline_threshold'] != 1e-8:
        raise ValueError('frozen five-track correction changed')
    counts = {'rows':0,'tested':0,'numerical_failures':0,'primary_legacy_p_below_2_5e_8':0,
              'v3_p_below_1e_8':0,'legacy_bh_q_below_0_05':0,'derived_five_track_q_below_0_05':0,
              'chromosome_order_violations':0,'adjacent_duplicate_keys':0}
    first = []
    last = []
    prev = None
    with gzip.open(ledger, 'rt', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != SOURCE_HEADER:
            raise ValueError('legacy 17-field schema changed')
        for row in reader:
            n = counts['rows'] + 1
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'malformed row {n}')
            if row['analysis_id'] != 'track-b-v1.0-pleiotropy' or row['pair_id'] != 'B':
                raise ValueError(f'wrong analysis/pair at row {n}')
            if row['analysis_status'] != 'TESTED' or row['numerical_error'] not in ('NA',''):
                raise ValueError(f'unsupported status/error mapping at row {n}')
            if int(row['within_pair_family_n']) != recovery['ledger_counts']['rows']:
                raise ValueError(f'legacy family denominator changed at row {n}')
            chrom, bp = int(row['CHR']), int(row['BP'])
            if not 1 <= chrom <= 22 or bp <= 0 or not row['SNP']:
                raise ValueError(f'invalid coordinate/SNP at row {n}')
            if row['A1'] not in 'ACGT' or row['A2'] not in 'ACGT' or row['A1'] == row['A2']:
                raise ValueError(f'invalid allele pair at row {n}')
            for field in ('Z1','Z2','T_PLACO_PLUS'):
                checked_float(row[field], field, n)
            for field in ('P1','P2','P_PLACO_PLUS','PLACO_BH_Q'):
                checked_float(row[field], field, n, probability=True)
            p = float(row['P_PLACO_PLUS'])
            q = float(row['PLACO_BH_Q'])
            key = (chrom,bp,row['SNP'])
            if prev is not None:
                counts['chromosome_order_violations'] += key < prev
                counts['adjacent_duplicate_keys'] += key == prev
            prev = key
            counts['rows'] = n
            counts['tested'] += 1
            counts['primary_legacy_p_below_2_5e_8'] += p <= 2.5e-8
            counts['v3_p_below_1e_8'] += p < 1e-8
            counts['legacy_bh_q_below_0_05'] += q < 0.05
            counts['derived_five_track_q_below_0_05'] += min(1.0, q*5) < 0.05
            mapped = {'SNP':row['SNP'],'CHR':row['CHR'],'BP':row['BP'],
                      'Z1':row['Z1'],'Z2':row['Z2'],'P_PLACO':row['P_PLACO_PLUS'],
                      'Q_WITHIN_PAIR':row['PLACO_BH_Q'],
                      'Q_BONFERRONI_ACROSS_PAIRS':repr(min(1.0,q*5)),
                      'status':'TESTED','headline':str(p < 1e-8)}
            if n <= 3:
                first.append(mapped)
            last.append(mapped)
            if len(last) > 3:
                last.pop(0)
    if counts['rows'] != recovery['ledger_counts']['rows'] or counts['primary_legacy_p_below_2_5e_8'] != recovery['ledger_counts']['primary']:
        raise ValueError('legacy rows/headline counts disagree with recovery audit')
    if counts['legacy_bh_q_below_0_05'] != recovery['ledger_counts']['bh']:
        raise ValueError('legacy BH count disagrees with recovery audit')
    if counts['chromosome_order_violations'] or counts['adjacent_duplicate_keys']:
        raise ValueError('ledger order/adjacent identity problem')
    output = {
        'analysis_id':'brain6-track-b-schema-crosswalk-v1',
        'scope':'READ_ONLY_TECHNICAL_CROSSWALK_NOT_CANONICAL_ADMISSION',
        'source_path':str(ledger),'source_sha256':expected_sha,
        'source_provenance_sha256':recovery['file_hashes']['archived_B_provenance']['sha256'],
        'target_template_path':str(TEMPLATE),'target_template_sha256':sha256(TEMPLATE),
        'family_lock_sha256':sha256(LOCK),'source_header':SOURCE_HEADER,'target_header':TARGET_HEADER,
        'field_map':{'SNP':'SNP','CHR':'CHR','BP':'BP','Z1':'Z1','Z2':'Z2',
                     'P_PLACO':'P_PLACO_PLUS','Q_WITHIN_PAIR':'PLACO_BH_Q',
                     'Q_BONFERRONI_ACROSS_PAIRS':'min(1, 5 * PLACO_BH_Q)',
                     'status':'analysis_status (TESTED only)',
                     'headline':'P_PLACO_PLUS < 1e-8'},
        'counts':counts,'first_three_proposed_rows':first,'last_three_proposed_rows':last,
        'admission_status':'BLOCKED_PENDING_EXPLICIT_FROZEN_RULE_AND_INDEPENDENT_VALIDATOR',
        'not_proven':['scientific equivalence of legacy and v3 implementation',
                      'frozen admission of V2 terminal-gate alias and V4 cleanup',
                      'independent whole-family correction/region eligibility']
    }
    output_path = OUT / 'schema_crosswalk_audit.json'
    stable_write(output_path, json.dumps(output, indent=2, sort_keys=True)+'\n')
    report = f'''# Protected Track B: technical legacy-to-v3 schema crosswalk

The archived B ledger is intact at SHA-256 `{expected_sha}`. A read-only scan checked all **{counts['rows']:,}** rows against its 17-field historical schema and the published v3 10-field template. Every row was `TESTED`, had the declared within-pair denominator, finite statistics and valid P/q values; the ledger was sorted with no adjacent duplicate key. This diagnostic did not create a candidate v3 output or import B.

| v3 field | Proposed historical source or deterministic expression |
|---|---|
'''
    for target, source in output['field_map'].items():
        report += f'| `{target}` | `{source}` |\n'
    report += f'''
The historical B headline threshold was P≤2.5e-8 across the original two primary Track B pairs ({counts['primary_legacy_p_below_2_5e_8']:,} rows). The frozen Brain6 v3 five-track headline is P<1e-8, which would select **{counts['v3_p_below_1e_8']:,}** B rows under the proposed adapter. The historical within-pair BH q<0.05 count is {counts['legacy_bh_q_below_0_05']:,}; multiplying that q by five and capping at one yields {counts['derived_five_track_q_below_0_05']:,} rows with proposed five-track q<0.05. These are **diagnostic counts, not admitted v3 family results**.

The conversion would copy SNP, chromosome, position and Z scores unchanged; rename `P_PLACO_PLUS` and `PLACO_BH_Q`; derive the five-track q and headline using the frozen v3 rules; and retain `TESTED`. The original A1/A2, P1/P2, statistic, error and family-size fields must remain preserved in the source ledger/provenance, even though the v3 display schema omits them. A read-only admission validator must additionally verify the exact historical gate/cleanup lineage, input/source identity, result hashes, and scientific compatibility before any import. No current frozen rule grants that admission; the protected slot remains blocked.
'''
    stable_write(OUT / 'schema_crosswalk_proposal.md', report)
    receipt = {'analysis_id':output['analysis_id'], 'source_sha256':expected_sha,
               'audit_sha256':sha256(output_path),'proposal_sha256':sha256(OUT/'schema_crosswalk_proposal.md'),
               'script_sha256':sha256(Path(__file__)), 'admission_status':output['admission_status']}
    stable_write(OUT / 'schema_crosswalk_receipt.json', json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'rows':counts['rows'],'v3_headline_diagnostic':counts['v3_p_below_1e_8'],
                      'admission_status':output['admission_status']}))

if __name__ == '__main__':
    main()
