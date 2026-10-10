#!/usr/bin/env python3
"""Additive table-routing and historical derived-column checks; no GWAS reads."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

P = Path(__file__).resolve().parents[1]
R = P/'reviews'
OUT = R/'independent_whole_extension_table_provenance_receipt_v4.json'
SEEN = {}


def text(path):
    path = Path(path)
    assert not path.name.endswith(('.sumstats.gz','.bgz','.bgz.partial'))
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
    info = {'sha256':h.hexdigest(),'bytes':path.stat().st_size}
    if str(path) in SEEN:
        assert SEEN[str(path)]==info
    SEEN[str(path)] = info
    return path.read_text()


def table(path):
    return list(csv.DictReader(text(path).splitlines(),delimiter='\t'))


def main():
    assert not OUT.exists()
    prior = json.loads(text(R/'independent_whole_extension_receipt_v4.json'))
    assert prior['status']=='PASS_WHOLE_EXTENSION_PROCESSED_INPUT_NATIVE_REPRODUCTION'
    pairs = table(P/'tables/extension_native_full_precision_rg_v4.tsv')
    h2 = table(P/'tables/extension_native_full_precision_h2_v4.tsv')
    indices = {}
    for row in pairs:
        path = Path(row['original_log'])
        if str(path) not in indices:
            lines = text(path).splitlines()
            start = lines.index('Summary of Genetic Correlation Results')
            index = {}
            for lineno in range(start+2,len(lines)):
                line = lines[lineno]
                if not line.strip():
                    break
                match = re.match(r'^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$',line)
                assert match
                key = Path(match[1]).name,Path(match[2]).name
                assert key not in index
                index[key] = lineno+1
            indices[str(path)] = index
        key = row['sleep_trait'],row['outcome_trait']
        assert row['stage']=='extension' and indices[str(path)][key]==int(row['original_log_line'])
    ratio_count = 0
    z_gate_changes,intercept_gate_changes = [],[]
    for row in h2:
        body = text(row['original_log'])
        assert row['stage']=='extension'
        source = re.search(r'^Read summary statistics for (\d+) SNPs\.',body,re.M)[1]
        regression = re.search(r'^After merging with regression SNP LD, (\d+) SNPs remain\.',body,re.M)[1]
        assert int(row['original_input_snp_count'])==int(source)
        assert int(row['original_ldsc_regression_snp_count'])==int(regression)
        oldz = float(row['original_h2'])/float(row['original_h2_se'])
        assert abs(oldz-float(row['original_h2_z']))<=max(1e-12*abs(oldz),1e-15)
        assert (row['original_pass_h2_z_ge_4']=='True')==(oldz>=4)
        assert (row['original_pass_intercept_le_1.2']=='True')==(float(row['original_LDSC_intercept'])<=1.2)
        ratio = re.search(r'^Ratio:\s*([\d.eE+\-]+)\s+\(([^)]+)\)',body,re.M)
        if ratio:
            ratio_count += 1
            assert row['original_attenuation_ratio_status']=='NUMERIC'
            assert format(float(row['original_attenuation_ratio']),'.4g')==format(float(ratio[1]),'.4g')
        else:
            assert 'Ratio: NA' in body or 'Ratio <' in body
        if row['h2_z_ge_4_full_precision']!=row['original_pass_h2_z_ge_4']:
            z_gate_changes.append(row['trait_id'])
        if row['intercept_le_1_2_full_precision']!=row['original_pass_intercept_le_1.2']:
            intercept_gate_changes.append(row['trait_id'])
    assert len(pairs)==1200 and len(h2)==100
    for path,info in list(SEEN.items()):
        text(path)
        assert SEEN[path]==info
    result = {'status':'PASS_ADDITIVE_EXTENSION_TABLE_PROVENANCE',
              'completed_utc':datetime.now(timezone.utc).isoformat(),
              'original_stock_log_line_routes_verified':len(pairs),'distinct_original_rg_logs':len(indices),
              'historical_h2_count_SNP_counts_derived_Z_and_flags_verified':len(h2),
              'historical_numeric_attenuation_ratios_verified':ratio_count,
              'native_vs_frozen_h2_Z_gate_changes':z_gate_changes,'native_vs_frozen_intercept_gate_changes':intercept_gate_changes,
              'GWAS_bodies_read':False,'native_workers_or_estimators':0,'consumed_artifacts_unchanged':SEEN,
              'scope':'Additive checks of routing/count/serialized derived columns; whole numerical/native provenance PASS remains in the separate hashed whole-extension receipt.'}
    with OUT.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['consumed_artifacts_unchanged','scope']}))


if __name__=='__main__':
    main()
