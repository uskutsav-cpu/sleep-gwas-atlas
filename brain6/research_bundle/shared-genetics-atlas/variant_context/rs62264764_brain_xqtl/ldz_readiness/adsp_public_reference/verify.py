#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'pair.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==1
x=r[0];assert x['#CHR']=='chr3' and x['POS1']=='117920728' and x['POS2']=='117923158' and x['REF:ALT1']=='G:A' and x['REF:ALT2']=='A:C'
a=float(x['R']);b=float(x['Rsq']);assert abs(a*a-b)<.00002
s=json.loads((p/'stream_receipt.json').read_text());assert s['compressed_bytes']==s['expected_bytes']==5043829695 and s['md5']==s['expected_md5']=='481f8ba5dcb3cd30c8c31bacc675a3cf'
manifest=list(csv.DictReader((p/'NG00067.v21.csv').open()));m=[v for v in manifest if v['file_name']=='adsp.wgs.r5.NHW.chr3.26042.2026.01.26.ld.output.txt.gz'];assert len(m)==1 and m[0]['md5']==s['md5'] and int(m[0]['file_size_bytes'])==s['compressed_bytes']
rec={'pair_rows':len(r),'R':a,'R2':b,'R_squared_from_rounded_R':a*a,'difference_from_reported_R2':abs(a*a-b),'manifest_match':True,'pair_sha256':hashlib.sha256((p/'pair.tsv').read_bytes()).hexdigest(),'note':'Integrity and arithmetic only; study-population representativeness, allele effect orientation and causality remain unverified.'}
(p/'verification.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
