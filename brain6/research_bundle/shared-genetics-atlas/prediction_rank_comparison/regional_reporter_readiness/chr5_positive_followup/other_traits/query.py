import csv,gzip,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;raw=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
sources={'mvp':(raw/'GCST90475826.tsv.gz','sha256','46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835',19703815),'hfrs':(raw/'GCST90707291.tsv','md5','38d9434f50a4493def9d2c9d1cc1b4b3',21294560)}
def digest(q,alg):
 h=hashlib.new(alg)
 with q.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
receipt={}
for label,(path,alg,expected,expected_n) in sources.items():
 h=digest(path,alg);assert h==expected;op=gzip.open if label=='mvp' else open;hits=[];n=0
 with op(path,'rt') as f:
  header=f.readline().rstrip('\r\n').split('\t');idx={k:header.index(k) for k in ('chromosome','base_pair_location','effect_allele','other_allele','p_value', 'rsid' if label=='mvp' else 'rs_id')};ik=idx['rsid' if label=='mvp' else 'rs_id']
  for line in f:
   n+=1;v=line.rstrip('\r\n').split('\t')
   if v[ik]=='rs12659965' or (v[idx['chromosome']]=='5' and v[idx['base_pair_location']]=='104568698'):hits.append(dict(zip(header,v)))
 assert n==expected_n
 exact=[r for r in hits if r['chromosome']=='5' and r['base_pair_location']=='104568698' and r['rsid' if label=='mvp' else 'rs_id']=='rs12659965' and {r['effect_allele'],r['other_allele']}=={'G','A'}]
 for x in exact:pv=float(x['p_value']);assert math.isfinite(pv) and 0<=pv<=1
 receipt[label]=dict(source=str(path),algorithm=alg,hash=h,rows=n,all_matches=hits,exact_matches=exact,status='MATCHED' if len(exact)==1 else ('NO_MATCH' if not exact else 'AMBIGUOUS'))
 print(label,'rows',n,'all matches',len(hits),'exact',len(exact),flush=True)
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
