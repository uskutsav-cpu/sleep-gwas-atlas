"""One streaming read of the authorized retained derivative; no raw/ref body."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import time

P = Path(__file__).resolve().parents[1]
plan = json.loads((P / 'manifests/finngen_preprocessing_plan_v4_3_8.json').read_text())
result_path = Path(plan['preprocessing_receipt'])
result = json.loads(result_path.read_text())
assert hashlib.sha256(result_path.read_bytes()).hexdigest() == '6d9a56ddf003c339ed7ccb50972b885d265d2265f6ede5b4266248610a5abe5b'
assert result['status'] == 'QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS'
assert result['source_before'] == result['source_after'] == plan['source_identity']
assert sum(result['terminal_dispositions'].values()) == result['source_rows'] == 21327010
assert result['derivative_rows'] == result['terminal_dispositions']['retained'] == 1156359
assert result['assumed_effective_N'] == plan['assumed_effective_N'] == 4 / (1 / 51643 + 1 / 446273)
assert result['diagnostic_counts']['supplied_P_finite_positive'] == result['derivative_rows']
assert result['diagnostic_counts']['P_Z_relative_difference_gt_0_1'] == 0
assert result['P_Z_diagnostics_used_for_filtering'] is False
assert result['scientific_source_admitted'] is result['independent_replication_established'] is False
assert result['verified_per_variant_N'] is result['verified_per_variant_INFO'] is False
assert result['full_source_gzip_CRC_and_EOF_verified'] is result['derivative_full_gzip_CRC_and_EOF_verified'] is True
path = Path(result['derivative'])
assert path == Path(plan['derivative']) and path.is_file() and not path.is_symlink()
assert not any(p.is_symlink() for p in path.parents)
started = time.monotonic()
compressed = hashlib.sha256()
decoded = hashlib.sha256()
identity = hashlib.sha256()
seen = set()  # Integer rsIDs only; bounded below128MiB measured RSS.
rows, compressed_bytes, decoded_bytes = 0, 0, 0
minimum_z, maximum_z = math.inf, -math.inf
class Reader:
    def __init__(self, handle): self.handle = handle
    def read(self, size=-1):
        global compressed_bytes
        value = self.handle.read(size)
        compressed.update(value)
        compressed_bytes += len(value)
        return value
    def tell(self): return self.handle.tell()
with path.open('rb') as raw:
    reader = Reader(raw)
    with gzip.GzipFile(fileobj=reader, mode='rb') as stream:
        header = stream.readline()
        assert header == b'SNP\tA1\tA2\tZ\tN\n'
        decoded.update(header); decoded_bytes += len(header)
        for line in stream:
            rows += 1
            decoded.update(line); decoded_bytes += len(line)
            fields = line.rstrip(b'\n').split(b'\t')
            assert len(fields) == 5 and line.endswith(b'\n') and b'\r' not in line
            snp, a1, a2, z_text, n_text = fields
            assert re.fullmatch(rb'rs[0-9]+', snp)
            number = int(snp[2:])
            assert number not in seen
            seen.add(number)
            assert a1 in [b'A', b'C', b'G', b'T'] and a2 in [b'A', b'C', b'G', b'T']
            assert a1 != a2 and frozenset([a1, a2]) not in [frozenset([b'A', b'T']), frozenset([b'C', b'G'])]
            z, n = float(z_text), float(n_text)
            assert math.isfinite(z) and math.isfinite(n) and n > 0
            assert n_text.decode() == format(plan['assumed_effective_N'], '.12g')
            minimum_z, maximum_z = min(minimum_z, z), max(maximum_z, z)
            identity.update(b'\t'.join(fields[:3]) + b'\n')
            if rows % 50000 == 0:
                assert time.monotonic() - started < 60
                assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < 128 << 20
    assert reader.read(1) == b''
assert compressed_bytes == path.stat().st_size
assert compressed.hexdigest() == result['derivative_sha256'] == 'e861abc30ece3ca4595e502082ca38d365093f0deaac558e5692fd97a9161092'
assert decoded.hexdigest() == result['derivative_decompressed_sha256']
assert identity.hexdigest() == result['ordered_SNP_allele_sha256']
assert rows == len(seen) == result['derivative_rows']
receipt = dict(schema='independent_actual_finngen_derivative_single_stream_v4_3_8', status='PASS',
    result_sha256=hashlib.sha256(result_path.read_bytes()).hexdigest(), derivative=str(path),
    derivative_sha256=compressed.hexdigest(), compressed_bytes=compressed_bytes,
    derivative_decompressed_sha256=decoded.hexdigest(), decompressed_bytes=decoded_bytes,
    ordered_SNP_allele_sha256=identity.hexdigest(), rows=rows, unique_rsids=len(seen),
    full_derivative_gzip_CRC_EOF_verified=True, all_Z_finite=True, all_N_literal_expected=True,
    all_alleles_single_distinct_nonambiguous=True, minimum_Z=minimum_z, maximum_Z=maximum_z,
    N_literal=format(plan['assumed_effective_N'], '.12g'), source_dispositions_sum=result['source_rows'],
    source_disposition_counts=result['terminal_dispositions'],
    source_P_and_coordinate_MAF_gates='Producer receipt and unchanged code/previous row-controls; not independently reread from raw.',
    raw_or_reference_body_reads=0, derivative_body_streams=1, real_workers_fits_network_or_locks=0,
    elapsed_seconds=time.monotonic()-started, max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
out = P / 'reviews/independent_finngen_actual_derivative_receipt_v4_3_8.json'
with out.open('x') as handle: json.dump(receipt, handle, indent=2); handle.write('\n')
print(json.dumps(receipt))
