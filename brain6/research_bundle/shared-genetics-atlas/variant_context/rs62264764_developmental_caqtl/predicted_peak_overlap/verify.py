import csv,gzip,hashlib,json,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[2].parent/'brain6-alphagenome-pilot/results/jbrowse_expanded';v='chr3:117920728:G>A';folder=root/'variants'/hashlib.sha256(v.encode()).hexdigest()[:20];m=json.loads((folder/'provenance.json').read_text());data=(folder/'track_window.f32.gz').read_bytes();assert hashlib.sha256(data).hexdigest()==m['track_window_sha256'];a=np.frombuffer(gzip.decompress(data),dtype='<f4').reshape(4,2,2001).astype('float64');i=[i for i,t in enumerate(m['tracks']) if t['output_type']=='dnase' and t['metadata'].get('biosample_life_stage')=='adult'][0];r=list(csv.DictReader((p/'segment_scores.tsv').open(),delimiter='\t'));origin=json.loads((p/'receipt.json').read_text())['track_window_start_0based'];maxerr=0
for x in r:
 start,end=int(x['start0'])-origin,int(x['end0'])-origin
 for j,key in [(0,'ref_sum'),(1,'alt_sum')]:maxerr=max(maxerr,abs(a[i,j,start:end].sum()-float(x[key])))
 assert end-start==int(x['bases'])
assert maxerr<1e-10
out={'numpy_track_sums_max_abs_error':maxerr,'rows':len(r),'segment_partition_matches':True,'meaning':'Numerical verification only; no evidence for measured genotype effects or matching tissue.'};(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
