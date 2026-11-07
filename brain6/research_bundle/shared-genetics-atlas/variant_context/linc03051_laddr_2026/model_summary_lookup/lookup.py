import csv, hashlib, io, json, re, zipfile
from pathlib import Path
import requests

HERE = Path(__file__).resolve().parent
URL = 'https://zenodo.org/api/records/19924423/files/twas_model_summaries.zip/content'
SIZE = 142716444

class Remote(io.RawIOBase):
    def __init__(self): self.pos = 0
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else SIZE + offset
        return self.pos
    def read(self, n=-1):
        if n < 0: n = SIZE - self.pos
        n = min(n, SIZE - self.pos)
        if not n: return b''
        r = requests.get(URL, headers={'Range': f'bytes={self.pos}-{self.pos+n-1}'}, params={'offset': self.pos}, timeout=60)
        r.raise_for_status()
        assert r.status_code == 206 and r.headers.get('Content-Range') == f'bytes {self.pos}-{self.pos+n-1}/{SIZE}' and len(r.content) == n
        self.pos += n
        return r.content

with (HERE.parent / 'adhd_unfiltered/matched_models.tsv').open() as f:
    targets = list(csv.DictReader(f, delimiter='\t'))
groups = {}
for row in targets:
    kind = {'latent_full':'ddp', 'latent_residual':'rddp'}[row['model']]
    groups.setdefault((kind, row['tissue']), []).append(row)
out, sources = [], []
with zipfile.ZipFile(Remote()) as archive:
    for (kind, tissue), rows in sorted(groups.items()):
        name = f'twas_model_summaries/{kind}-{tissue}.profile'
        info = archive.getinfo(name)
        data = archive.read(name)  # zipfile checks CRC on complete read
        sources.append({'member':name, 'uncompressed_bytes':len(data), 'zip_crc32':f'{info.CRC:08x}', 'sha256':hashlib.sha256(data).hexdigest()})
        records = list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))
        for row in rows:
            hits = [x for x in records if x['id'] == row['feature']]
            assert len(hits) == 1, (name, row['feature'], len(hits))
            h = {k: v.strip() for k,v in hits[0].items()}
            method = row['weights_model']
            assert re.fullmatch(r'top1|blup|enet|bslmm|lasso', method)
            out.append({'tissue':tissue,'model':row['model'],'feature':row['feature'],'selected_method':method,
                        'reported_nsnps':h['nsnps'],'reported_hsq':h['hsq'],'reported_hsq_se':h['hsq.se'],
                        'reported_method_r2':h[method+'.r2'],'reported_method_p':h[method+'.pv'],
                        'adhd_p':row['adhd_p'],'sleep_p':row['sleep_p']})
assert len(out) == len(targets) == 8
with (HERE/'selected_profiles.tsv').open('w') as f:
    w=csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)
(HERE/'receipt.json').write_text(json.dumps({'archive':URL,'archive_bytes':SIZE,'publisher_archive_md5':'47fd6eb2192a891a3946c9bcda20963e','full_archive_verified':False,'members':sources,'selected_rows':len(out),'limits':'ZIP CRC checks transport consistency, not full archive authenticity; profiles are not model weights.'},indent=2)+'\n')
print(f'complete: {len(out)} selected profiles; {len(sources)} archive members')
