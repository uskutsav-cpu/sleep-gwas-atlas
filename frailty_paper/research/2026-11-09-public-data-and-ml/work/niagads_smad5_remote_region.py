#!/usr/bin/env python3
"""Extract a tabix-indexed region from an uncompressed remote tar via byte ranges."""
import csv, gzip, hashlib, io, json, shutil, struct, subprocess, tempfile, sys
from pathlib import Path
import pysam

root = Path(__file__).resolve().parent.parent
manifest = root / 'work/niagads_snuc_all_manifest.txt'
url = 'https://st1.niagads.org/portal/v1/download-public/NG00184/fileset/57103'
assert len(sys.argv) <= 2 and (len(sys.argv) == 1 or sys.argv[1] == '--tad')
is_tad = len(sys.argv) == 2
start, stop = (134640000, 137760000) if is_tad else (135132844, 137132845)  # zero-based half-open GRCh38
suffix = '_tad' if is_tad else ''
offset = 0
data = index = None
with manifest.open() as f:
    for row in csv.DictReader(f, delimiter='\t'):
        name, size = row['File name'], int(row['File size'])
        if 'ROSMAP_CUIMC1_2_MIT.snuc-eQTL.exc.' in name and '_v1_5_snp_all.' in name:
            entry = {'name': name, 'size': size, 'md5': row['md5Sum'], 'header_offset': offset, 'data_offset': offset + 1536}
            if name.endswith('.bed.gz'): data = entry
            if name.endswith('.bed.gz.tbi'): index = entry
        # These GNU tar members have 512-byte headers + 512-byte long names + a 512-byte file header.
        offset += 1536 + (size + 511) // 512 * 512
assert data and index and index['data_offset'] > data['data_offset']

def fetch(lo, hi):
    proc = subprocess.run(['curl', '-fsSL', '--max-time', '40', '-r', f'{lo}-{hi}', url],
                          capture_output=True, check=True, timeout=45)
    if len(proc.stdout) != hi - lo + 1:
        raise ValueError(f'Incorrect byte range size for {lo}-{hi}: {len(proc.stdout)}')
    return proc.stdout

# Verify expected long-name header at the calculated location, then verify complete index against manifest.
header = fetch(data['header_offset'], data['data_offset'] - 1)
assert header[:13] == b'././@LongLink' and data['name'].encode() in header
idx = fetch(index['data_offset'], index['data_offset'] + index['size'] - 1)
assert hashlib.md5(idx).hexdigest() == index['md5']
tb = gzip.decompress(idx)
assert tb[:4] == b'TBI\1'
n_ref, fmt, seq, cbeg, cend, meta, skip, name_len = struct.unpack_from('<8i', tb, 4)
assert n_ref == 1 and tb[36:36 + name_len] == b'chr5\0'
p = 36 + name_len
n_bins = struct.unpack_from('<i', tb, p)[0]; p += 4
bins = {0}
for shift, base in [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)]:
    bins.update(range(base + (start >> shift), base + ((stop - 1) >> shift) + 1))
ranges = [(0, min(data['size'] - 1, 65535))]
for _ in range(n_bins):
    bin_id, n_chunks = struct.unpack_from('<Ii', tb, p); p += 8
    for _ in range(n_chunks):
        a, b = struct.unpack_from('<QQ', tb, p); p += 16
        if bin_id in bins:
            ranges.append((a >> 16, min(data['size'] - 1, (b >> 16) + 65535)))
ranges.sort(); merged = []
for lo, hi in ranges:
    if merged and lo <= merged[-1][1] + 1:
        merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
    else:
        merged.append((lo, hi))
assert len(merged) > 1 and sum(hi - lo + 1 for lo, hi in merged) < 100_000_000
with tempfile.TemporaryDirectory(dir=root / 'work') as tmp:
    file = Path(tmp) / 'subset.bed.gz'
    with file.open('wb') as f:
        f.truncate(data['size'])  # Sparse file: only fetched blocks occupy disk space.
        for lo, hi in merged:
            f.seek(lo); f.write(fetch(data['data_offset'] + lo, data['data_offset'] + hi))
    Path(str(file) + '.tbi').write_bytes(idx)
    with pysam.TabixFile(str(file)) as tabix:
        columns = list(tabix.header)[-1].lstrip('#').split('\t')
        lines = list(tabix.fetch('chr5', start, stop))
    rows = [dict(zip(columns, line.split('\t'))) for line in lines]
    filtered = [r for r in rows if r.get('target_ensembl_id') == 'ENSG00000113658']
    assert filtered, 'No SMAD5 rows; check schema and query'
    assert len({r['variant_id'] for r in filtered}) == len(filtered), 'Duplicate target variant rows'
    output = root / f'outputs/stage6_niagads_smad5_exc_region{suffix}.tsv.gz'
    with gzip.open(output, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t'); writer.writeheader(); writer.writerows(filtered)
receipt = {'source': url, 'archive_member': data, 'index_member': index, 'index_md5_verified': True,
           'archive_header_verified': True, 'archive_full_md5_verified': False,
           'GRCh38_region_zero_based_half_open': [start, stop], 'compressed_byte_ranges': merged,
           'compressed_bytes_downloaded': sum(hi-lo+1 for lo,hi in merged),
           'all_target_rows_in_window': len(rows), 'smad5_rows': len(filtered),
           'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
           'limitations': 'Indexed regional access only; no full member checksum; target-specific model and phenotype harmonization needed.'}
(root / f'outputs/stage6_niagads_smad5_exc_region{suffix}_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
