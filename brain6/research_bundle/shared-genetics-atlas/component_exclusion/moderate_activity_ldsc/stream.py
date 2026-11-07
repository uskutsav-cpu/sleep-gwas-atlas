import gzip, hashlib, json, os, sys, urllib.request
from pathlib import Path
url='https://www.dropbox.com/s/foxe08tsyxi7bzq/884.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1'
with urllib.request.urlopen(url, timeout=120) as r: data=r.read()
plain=gzip.decompress(data)
lines=plain.splitlines()
header=lines[0].decode().split('\t')
if header != ['SNP','A1','A2','N','Z']: raise ValueError(f'Unexpected header: {header}')
if len(lines)<1000000: raise ValueError('Unexpectedly few rows')
Path(__file__).with_name('receipt_'+os.environ['STREAM_TAG']+'.json').write_text(json.dumps({'url':url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'rows':len(lines)-1,'header':header},indent=2)+'\n')
sys.stdout.buffer.write(plain)
