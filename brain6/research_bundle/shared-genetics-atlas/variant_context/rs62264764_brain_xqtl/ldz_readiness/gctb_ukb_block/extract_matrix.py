#!/usr/bin/env python3
"""Extract the allele-aligned matched submatrix, checking the whole ZIP member CRC."""
import collections, csv, hashlib, json, zipfile
from pathlib import Path
import numpy as np
from remote_zip import Remote

p = Path(__file__).resolve().parent
meta = json.loads((p / 'archive_inventory.json').read_text())
source = p.parents[2] / 'rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
with source.open() as f:
    variants = list(csv.DictReader(f, delimiter='\t'))
with (p / 'block341_snp.tsv').open() as f:
    snps = list(csv.DictReader(f, delimiter='\t'))
assert len(snps) == 13280
bypos = collections.defaultdict(list)
for i, s in enumerate(snps):
    bypos[int(s['PhysPos'])].append((i, s))
selected = []
for v in variants:
    hits = [(i, s) for i, s in bypos[int(v['pos37'])]
            if {s['A1'], s['A2']} == {v['ref'], v['alt']}]
    if len(hits) == 1:
        i, s = hits[0]
        selected.append((i, v, s, 1 if s['A1'] == v['alt'] else -1))
assert len({x[0] for x in selected}) == len(selected)
selected.sort(key=lambda x: x[0])
assert len(selected) > 1000
indices = np.array([x[0] for x in selected])
signs = np.array([x[3] for x in selected], dtype=np.float32)
result = np.empty((len(selected), len(selected)), dtype=np.float32)
remote = Remote(meta['url'], meta['archive_bytes'])
with zipfile.ZipFile(remote) as archive:
    info = archive.getinfo('ldm13M/block341.ldm.bin')
    assert info.file_size == 4 * len(snps) * (len(snps) + 1) // 2
    with archive.open(info) as member:
        k = 0
        for row in range(len(snps)):
            data = member.read(4 * (row + 1))
            if len(data) != 4 * (row + 1): raise IOError('truncated matrix')
            if k < len(indices) and indices[k] == row:
                values = np.frombuffer(data, dtype='<f4')[indices[:k+1]] * signs[k] * signs[:k+1]
                result[k, :k+1] = values
                result[:k+1, k] = values
                k += 1
        assert member.read(1) == b'' and k == len(indices)
assert np.isfinite(result).all()
np.savez_compressed(p / 'matched_ld.npz', ld=result, variant=np.array([x[1]['variant'] for x in selected]))
with (p / 'matched_snps.tsv').open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['variant','pos37','ref','alt','gctb_id','gctb_a1','gctb_a2','gctb_a1freq','row_index','alignment_sign'])
    for i, v, s, sign in selected:
        w.writerow([v['variant'],v['pos37'],v['ref'],v['alt'],s['ID'],s['A1'],s['A2'],s['A1Freq'],i,sign])
receipt = {'source_url':meta['url'],'member':info.filename,'zip_member_crc32':f'{info.CRC:08x}','zip_member_crc_passed':True,'member_uncompressed_bytes':info.file_size,'http_range_bytes_transferred':remote.transferred,'source_variants':len(variants),'matched_variants':len(selected),'missing_or_ambiguous':len(variants)-len(selected),'matrix_sha256':hashlib.sha256((p/'matched_ld.npz').read_bytes()).hexdigest(),'matrix_diagonal_min':float(np.diag(result).min()),'matrix_diagonal_max':float(np.diag(result).max()),'matrix_min':float(result.min()),'matrix_max':float(result.max()),'limitations':'ZIP member CRC checked; complete archive not downloaded or hashed. Match by GRCh37 position and exact allele set; analysis SNPs missing from the reference were excluded.'}
(p/'matrix_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
