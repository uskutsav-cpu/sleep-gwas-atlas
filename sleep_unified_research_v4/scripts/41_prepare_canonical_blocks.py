#!/usr/bin/env python3
"""Freeze reference-only genomic intervals; no GWAS/estimator outcome is read."""
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import resource
import shutil

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def guards():
    assert shutil.disk_usage('/System/Volumes/Data').free >= 3 * 1024**3
    assert shutil.disk_usage(SSD).free >= 5 * 1024**3
    # macOS reports ru_maxrss in bytes.
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < 768 * 1024**2


def main():
    plan_path = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
    plan = json.loads(plan_path.read_text())
    dependencies = plan['dependencies_sha256']
    refs = {str(p): h for p, h in dependencies.items() if p.endswith('.l2.ldscore.gz')}
    assert len(refs) == 22
    coordinate = json.loads((PACKAGE / 'statistical_validation/genomicsem_specialist_forensic_checks_v4.json').read_text())['existing_coordinate_map']
    coord_path = Path(coordinate['path'])
    coord_hash = '6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4'
    assert sha(coord_path) == coord_hash
    destination = PACKAGE / 'statistical_validation/canonical_200_intervals_v4.tsv'
    receipt_path = PACKAGE / 'manifests/canonical_200_interval_preparation_v4.json'
    assert not destination.exists() and not receipt_path.exists()
    guards()
    chromosome_positions, seen = {}, set()
    all_ordered, construction_ordered = hashlib.sha256(), hashlib.sha256()
    excluded, counts_all, reference_headers = Counter(), {}, {}
    for chromosome in range(1, 23):
        path = next(Path(p) for p in refs if Path(p).name == f'{chromosome}.l2.ldscore.gz')
        assert sha(path) == refs[str(path)]
        positions, previous, count = [], -1, 0
        with gzip.open(path, 'rt') as f:
            header = f.readline().split()
            assert {'CHR', 'SNP', 'BP'} <= set(header)
            ci, si, bi = [header.index(k) for k in ('CHR', 'SNP', 'BP')]
            reference_headers[str(path)] = header
            for line in f:
                vals = line.split(); c, snp, bp = int(vals[ci]), vals[si], int(vals[bi])
                assert c == chromosome and bp >= previous and bp >= 1
                assert snp not in seen
                seen.add(snp); previous = bp; count += 1
                encoded = f'{c}\t{snp}\t{bp}\n'.encode()
                all_ordered.update(encoded)
                if c == 6 and 25_000_000 <= bp <= 34_000_000:
                    excluded['original_extended_MHC'] += 1
                    continue
                positions.append(bp); construction_ordered.update(encoded)
        chromosome_positions[chromosome] = positions
        counts_all[chromosome] = count
        guards()
    counts = {c: len(ps) for c, ps in chromosome_positions.items()}
    total = sum(counts.values())
    assert all(counts.values())
    # Hamilton apportionment: one guaranteed block per chromosome, 178 by size.
    quotient = {c: 178 * counts[c] // total for c in counts}
    remainder = {c: 178 * counts[c] % total for c in counts}
    allocation = {c: 1 + quotient[c] for c in counts}
    outstanding = 200 - sum(allocation.values())
    for c in sorted(counts, key=lambda c: (-remainder[c], c))[:outstanding]:
        allocation[c] += 1
    assert sum(allocation.values()) == 200
    intervals, block_id, construction_counts = [], 0, []
    for c in range(1, 23):
        ps, n, b = chromosome_positions[c], counts[c], allocation[c]
        cuts = [1]
        for j in range(1, b):
            target = j * n // b
            bp = ps[target]
            assert bp > cuts[-1], 'TIED_QUANTILE_COORDINATE_REQUIRES_RESULT_FREE_AMENDMENT'
            cuts.append(bp)
        cuts.append(max(ps) + 1)
        observed = [0] * b
        which = 0
        for bp in ps:
            while bp >= cuts[which + 1]:
                which += 1
            observed[which] += 1
        assert all(observed)
        for j in range(b):
            intervals.append({'block_id': block_id, 'chromosome': c, 'start_inclusive': cuts[j],
                              'end_exclusive': cuts[j + 1], 'construction_reference_snps': observed[j]})
            construction_counts.append(observed[j]); block_id += 1
    assert len(intervals) == 200
    with destination.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(intervals[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(intervals)
    assert all(sha(p) == h for p, h in refs.items()) and sha(coord_path) == coord_hash
    guards()
    receipt = {'recorded_utc': datetime.now(timezone.utc).isoformat(), 'build': 'GRCh37/hg19', 'block_count': 200,
               'interval_tsv_path': str(destination), 'interval_tsv_sha256': sha(destination),
               'reference_sha256': refs, 'coordinate_map_path': str(coord_path), 'coordinate_map_sha256': coord_hash,
               'native_plan_sha256': sha(plan_path), 'builder_sha256': sha(Path(__file__)),
               'ordered_all_reference_records_sha256': all_ordered.hexdigest(),
               'ordered_construction_records_sha256': construction_ordered.hexdigest(),
               'all_reference_snp_count': sum(counts_all.values()), 'construction_reference_snp_count': total,
               'chromosome_counts_all': counts_all, 'chromosome_counts_construction': counts,
               'chromosome_block_allocation': allocation, 'block_reference_snp_counts': construction_counts,
               'reference_headers': reference_headers, 'static_construction_exclusions': dict(excluded),
               'partition_amendment': 'Outcome-blind original extended MHC exclusion applied to reference construction counts to avoid intervals populated solely by scientifically excluded MHC SNPs; no GWAS retained sets used.',
               'allocation_rule': 'One block per chromosome plus Hamilton allocation of178 by reference SNP count; fractional ties numeric chromosome.',
               'cut_rule': 'For within-chromosome block j use BP at floor(j*N/B); equal BP values never split. Start1, end max construction BP+1. Reject repeated cut BP.',
               'all_source_reference_hashes_unchanged_after': True, 'GWAS_outcomes_read': False,
               'estimator_calls': 0, 'final_pair_coordinate_coverage_verified': False,
               'empirical_uncertainty_calibrated': False, 'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    folder = SSD / 'canonical_calibration'; folder.mkdir(exist_ok=True)
    (folder / destination.name).write_bytes(destination.read_bytes())
    (folder / receipt_path.name).write_bytes(receipt_path.read_bytes())
    print(json.dumps({k: receipt[k] for k in ('block_count', 'interval_tsv_sha256', 'all_reference_snp_count', 'construction_reference_snp_count', 'peak_rss_bytes', 'estimator_calls')}))


if __name__ == '__main__':
    main()
