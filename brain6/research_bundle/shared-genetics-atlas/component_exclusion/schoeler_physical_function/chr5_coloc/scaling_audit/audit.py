import csv,hashlib,json,math,statistics
from pathlib import Path
p=Path(__file__).resolve().parent;region=p.parent
path=region/'insomnia_matched.tsv';receipt=json.loads((region/'match_receipt.json').read_text())
hashv=hashlib.sha256(path.read_bytes()).hexdigest();assert hashv==receipt['results']['insomnia']['matched_sha256']
with (region/'common_support/positions.tsv').open() as f:positions={int(x['pos38']) for x in csv.DictReader(f,delimiter='\t')}
with path.open() as f:rows=[x for x in csv.DictReader(f,delimiter='\t') if int(x['pos38']) in positions]
assert len(rows)==len(positions)==2359
D=[];S=[]
for x in rows:
 b,se,pval,freq=map(float,(x['physical_beta'],x['physical_se'],x['physical_p'],x['physical_frequency']))
 assert all(map(math.isfinite,(b,se,pval,freq))) and se>0 and 0<pval<=1 and 0<freq<1
 pred=math.erfc(abs(b/se)/math.sqrt(2));assert pred>0
 D.append(abs(math.log10(pval)-math.log10(pred)));S.append(se*math.sqrt(2*405979*freq*(1-freq)))
out={'input_sha256':hashv,'n_markers':len(rows),'N_assumed_for_proxy':405979,'abs_log10_p_difference':{'median':statistics.median(D),'max':max(D)},'simplified_scale_proxy':{'median':statistics.median(S),'min':min(S),'max':max(S)},'warning':'Proxy is conditional on simplified linear-regression assumptions, not a measured or validated trait SD; do not tune prior from this audit.'}
(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
