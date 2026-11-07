#!/usr/bin/env python3
import collections,csv,hashlib,json,math
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[3]
prior=json.loads((p.parent/'matching_receipt.json').read_text())['pairs']['exc_insomnia'];qtl=Path(prior['file'])
meta=json.loads((atlas/'source_provenance/raw_source_hashes.json').read_text())['FinnGen_R13_insomnia_GRCh38'];source=Path(meta['path']);idx=p/'finngen.tbi'
def digest(x,algo='sha256'):
 h=hashlib.new(algo)
 with x.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
assert digest(qtl)==prior['sha256'] and digest(source)==meta['sha256'] and digest(idx,'md5')=='e275bcbdd87ff79f41ec03294e44c1b1'
with qtl.open() as f:q=list(csv.DictReader(f,delimiter='\t'))
lookup=collections.defaultdict(list)
with pysam.TabixFile(str(source),index=str(idx)) as tb:
 for line in tb.fetch('3',117423157,118423158):
  x=line.split('\t');assert len(x)==13
  if len(x[2])==len(x[3])==1 and x[2] in 'ACGT' and x[3] in 'ACGT':lookup[(int(x[1]),frozenset((x[2],x[3])))].append(x)
rows=[];counts=collections.Counter()
for x in q:
 key=(int(x['pos38']),frozenset((x['ref'],x['alt'])));found=lookup.get(key,[])
 if len(found)!=1:counts['no_unique_finngen']+=1;continue
 y=found[0]
 try:
  beta,se,pv,af=map(float,(y[8],y[9],y[6],y[10]));qaf=float(x['qtl_alt_freq'])
  assert all(map(math.isfinite,(beta,se,pv,af,qaf))) and se>0 and 0<=pv<=1 and 0<af<1
  if y[3]!=x['alt']:beta=-beta;af=1-af
  rows.append({**x,'fg_rsid':y[4],'fg_beta_alt':beta,'fg_se':se,'fg_p':pv,'fg_alt_freq':af,'abs_af_difference':abs(qaf-af)});counts['matched']+=1
 except (ValueError,AssertionError):counts['invalid']+=1
assert len(rows)>=500 and len({x['variant'] for x in rows})==len(rows)
with (p/'matched.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
a=sorted(float(x['abs_af_difference']) for x in rows)
rec={'qtl_matched_input':str(qtl),'qtl_sha256':digest(qtl),'finngen_source':str(source),'finngen_sha256':digest(source),'index_md5':digest(idx,'md5'),'index_sha256':digest(idx),'index_generation':'1777989524607662','input_count':len(q),'counts':dict(counts),'median_absolute_af_difference':a[len(a)//2],'p95_absolute_af_difference':a[int(.95*(len(a)-1))],'maximum_absolute_af_difference':max(a),'input_source_note':'FinnGen is R13 registry insomnia; QTL source is a selected partial excitatory brain extract.'}
(p/'matching_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
