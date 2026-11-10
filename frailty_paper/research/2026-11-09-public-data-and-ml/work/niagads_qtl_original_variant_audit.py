import csv,gzip,re,json,hashlib
from collections import Counter
from pathlib import Path
q=Path('outputs/stage6_niagads_smad5_exc_region.tsv.gz'); a=Path('outputs/stage6_niagads_smad5_exc_Mak_FI_aligned.csv')
ids={r['variant_GRCh38'] for r in csv.DictReader(a.open())}; rows=[]
with gzip.open(q,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['variant_id'] not in ids: continue
  user=dict(x.split('=',1) for x in r['user_input'].split(';') if '=' in x)
  orig=user.get('original_variant_id',''); m=re.fullmatch(r'chr(\w+):(\d+)_([ACGT]+)_([ACGT]+)',orig)
  if not m: status='unparsed'
  else:
   ch,pos,ref,alt=m.groups(); p=int(pos)
   while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
   while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:ref=ref[1:];alt=alt[1:];p+=1
   status='normalized_matches' if f'chr{ch}:{p}:{ref}:{alt}'==r['variant_id'] else 'normalized_differs'
  rows.append({'variant':r['variant_id'],'original':orig,'status':status,'reported_rsid':dict(x.split('=',1) for x in r['QC_info'].split(';') if '=' in x).get('dbSNP_rsid')})
assert len(rows)==len(ids)
c=Counter(x['status'] for x in rows); mult=[x for x in rows if len(x['original'].split('_')[-1])>1 or len(x['original'].split('_')[-2])>1]
out={'source_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'alignment_sha256':hashlib.sha256(a.read_bytes()).hexdigest(),'n':len(rows),'statuses':dict(c),'multibase_original_records':len(mult),'nonmatching_records':[x for x in rows if x['status']!='normalized_matches'],'multibase_examples':mult[:8],'limitations':'String-level minimal common-prefix/suffix normalization only; does not authenticate reference sequence, phase, source AF or cross-study alleles.'}
Path('outputs/stage6_niagads_qtl_original_variant_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='nonmatching_records'},indent=2));print('mismatch examples',out['nonmatching_records'][:8])
