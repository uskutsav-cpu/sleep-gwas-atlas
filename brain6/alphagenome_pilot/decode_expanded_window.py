#!/usr/bin/env python3
"""Decode one saved variant window into a human-readable TSV."""
import argparse,csv,gzip,hashlib,json,struct
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('variant_folder',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
m=json.loads((a.variant_folder/'provenance.json').read_text());b=(a.variant_folder/'track_window.f32.gz').read_bytes()
if hashlib.sha256(b).hexdigest()!=m['track_window_sha256']:raise ValueError('Compressed output hash mismatch')
v=struct.unpack('<16008f',gzip.decompress(b));pos=m['variant']['position'];a.out.parent.mkdir(parents=True,exist_ok=True)
with a.out.open('x',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['chromosome','position_grch38_1based','offset_from_variant_bp','output_type','track_ordinal','name','life_stage','data_source','strand','reference','alternate','delta_alt_minus_ref'])
 for i,t in enumerate(m['tracks']):
  d=t['metadata'];r=v[i*4002:i*4002+2001];alt=v[i*4002+2001:(i+1)*4002]
  for j,(x,y) in enumerate(zip(r,alt)):w.writerow([m['variant']['chromosome'],pos+j-1000,j-1000,t['output_type'],t['track_ordinal'],d.get('name'),d.get('biosample_life_stage'),d.get('data_source'),d.get('strand'),repr(x),repr(y),repr(y-x)])
print('Wrote 8,004 track rows')
