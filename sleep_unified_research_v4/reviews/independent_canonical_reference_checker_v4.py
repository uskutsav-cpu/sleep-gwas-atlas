#!/usr/bin/env python3
"""Outcome-free streaming validation of reference-only canonical intervals."""
from collections import Counter
import csv
import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import resource
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
SEEN = {}
STARTED = time.monotonic()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
    return h.hexdigest()


def record(path):
    SEEN[str(path)] = {'sha256':sha(path),'bytes':path.stat().st_size}
    return path


def guard():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        rss *= 1024
    assert rss < 128*1024**2 and shutil.disk_usage(ROOT).free >= 128*1024**2
    assert time.monotonic()-STARTED < 600
    return rss


def main():
    assert shutil.disk_usage(ROOT).free >= 256*1024**2
    record(Path(__file__))
    builder = record(P/'scripts/41_prepare_canonical_blocks.py')
    receipt = json.loads(record(P/'manifests/canonical_200_interval_preparation_v4.json').read_text())
    intervals = list(csv.DictReader(record(Path(receipt['interval_tsv_path'])).open(),delimiter='\t'))
    intervals = [{k:int(v) for k,v in r.items()} for r in intervals]
    assert sha(builder) == receipt['builder_sha256']
    assert sha(Path(receipt['interval_tsv_path'])) == receipt['interval_tsv_sha256']
    assert [r['block_id'] for r in intervals] == list(range(200))
    assert receipt['estimator_calls'] == 0 and receipt['GWAS_outcomes_read'] is False
    assert receipt['final_pair_coordinate_coverage_verified'] is False and receipt['empirical_uncertainty_calibrated'] is False
    assert sha(record(Path(receipt['coordinate_map_path']))) == receipt['coordinate_map_sha256']
    nplan = P/'manifests/ssd_native_execution_plan_v4_3.json'
    assert sha(record(nplan)) == receipt['native_plan_sha256']
    nd = json.loads(nplan.read_text())
    assert receipt['reference_sha256'] == {s:h for s,h in nd['dependencies_sha256'].items() if s.endswith('.l2.ldscore.gz')}
    all_stream,construction_stream = hashlib.sha256(),hashlib.sha256()
    counts,all_counts,allocations,block_counts,exclusions,checks = {},{}, {}, {},Counter(),[]
    for chromosome in range(1,23):
        guard()
        path = next(Path(s) for s in receipt['reference_sha256'] if Path(s).name == f'{chromosome}.l2.ldscore.gz')
        assert sha(record(path)) == receipt['reference_sha256'][str(path)]
        blocks = [r for r in intervals if r['chromosome'] == chromosome]
        assert blocks and blocks[0]['start_inclusive'] == 1
        assert all(a['end_exclusive'] == b['start_inclusive'] for a,b in zip(blocks,blocks[1:]))
        assert all(r['start_inclusive'] < r['end_exclusive'] for r in blocks)
        positions = [];observed = [0]*len(blocks);which = 0;prior = 0;total = 0
        with gzip.open(path,'rt') as f:
            header = f.readline().split();cols = {k:header.index(k) for k in ('CHR','SNP','BP')}
            for line in f:
                row = line.split();c,snp,bp = int(row[cols['CHR']]),row[cols['SNP']],int(row[cols['BP']])
                assert c == chromosome and bp >= prior and bp >= 1
                prior = bp;total += 1
                encoded = f'{c}\t{snp}\t{bp}\n'.encode();all_stream.update(encoded)
                if c == 6 and 25000000 <= bp <= 34000000:
                    exclusions['original_extended_MHC'] += 1
                    continue
                construction_stream.update(encoded);positions.append(bp)
                while which+1 < len(blocks) and bp >= blocks[which]['end_exclusive']:
                    which += 1
                assert blocks[which]['start_inclusive'] <= bp < blocks[which]['end_exclusive']
                observed[which] += 1
        n = len(positions);b = len(blocks)
        assert blocks[-1]['end_exclusive'] == max(positions)+1
        expected_cuts = [1]+[positions[(j*n)//b] for j in range(1,b)]+[max(positions)+1]
        assert [r['start_inclusive'] for r in blocks]+[blocks[-1]['end_exclusive']] == expected_cuts
        assert observed == [r['construction_reference_snps'] for r in blocks] and all(observed)
        all_counts[str(chromosome)] = total;counts[str(chromosome)] = n;allocations[str(chromosome)] = b
        block_counts.update({str(r['block_id']):v for r,v in zip(blocks,observed)})
        checks.append({'chromosome':chromosome,'all_reference_rows':total,'construction_rows':n,'blocks':b,'all_rows_in_exactly_one_interval':True,'quantile_cuts_reproduced':True,'membership_counts_reproduced':True})
    total = sum(counts.values())
    # Independently use floating Hamilton quotas; these small integer products
    # are exactly representable and no quota lies at an integer/tie boundary.
    quota = {c:178*counts[c]/total for c in counts}
    expected = {c:1+math.floor(q) for c,q in quota.items()}
    remaining = 200-sum(expected.values())
    order = sorted(counts,key=lambda c: (-(quota[c]-math.floor(quota[c])),int(c)))
    for c in order[:remaining]:
        expected[c] += 1
    assert expected == allocations == receipt['chromosome_block_allocation']
    assert counts == receipt['chromosome_counts_construction'] and all_counts == receipt['chromosome_counts_all']
    assert sum(all_counts.values()) == receipt['all_reference_snp_count'] and total == receipt['construction_reference_snp_count']
    assert [block_counts[str(i)] for i in range(200)] == receipt['block_reference_snp_counts']
    assert dict(exclusions) == receipt['static_construction_exclusions']
    assert all_stream.hexdigest() == receipt['ordered_all_reference_records_sha256']
    assert construction_stream.hexdigest() == receipt['ordered_construction_records_sha256']
    unchanged = all(sha(Path(s)) == v['sha256'] for s,v in SEEN.items());assert unchanged
    out = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'REFERENCE_ONLY_ALLOCATION_AND_MEMBERSHIP_PASS',
           'build':receipt['build'],'block_count':200,'reference_rows':sum(all_counts.values()),'construction_rows':total,
           'reference_MHC_exclusion_count':sum(exclusions.values()),'reference_already_lacks_extended_MHC':not exclusions,
           'allocation':expected,'checks':checks,'ordered_all_reference_sha256':all_stream.hexdigest(),'ordered_construction_sha256':construction_stream.hexdigest(),
           'all_rows_in_exactly_one_interval':True,'all_source_hashes_unchanged':unchanged,'inputs':SEEN,
           'coordinate_map_content_alignment_independently_checked':False,'global_reference_rsID_uniqueness_independently_checked':False,
           'genome_length_coverage_claimed':False,'final_pair_coordinate_coverage_verified':False,'cross_fit_delete_array_alignment_certified':False,
           'empirical_uncertainty_calibrated':False,'GWAS_source_or_estimator_outcomes_read':False,'native_estimator_called':False,
           'resource_limits':{'rss_cap_bytes':128*1024**2,'hash_buffer_bytes':65536,'internal_launch_bytes':256*1024**2,'internal_emergency_bytes':128*1024**2,'seconds_limit':600},
           'resource_usage':{'peak_reviewer_rss_bytes':guard(),'elapsed_seconds':time.monotonic()-STARTED}}
    dest = P/'reviews/independent_canonical_reference_receipt_v4.json';assert not dest.exists()
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['inputs','checks']},indent=2))


if __name__ == '__main__':
    main()
