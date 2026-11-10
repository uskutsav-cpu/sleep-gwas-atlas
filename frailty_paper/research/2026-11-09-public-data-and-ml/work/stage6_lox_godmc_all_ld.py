import json,csv
import numpy as np
from pathlib import Path
rows=json.load(open('work/godmc_cpg_cg12594244.json'))['assoc_meta'];lookup={}
for r in rows:
 try:
  _,chrom,pos,kind=r['snp'].replace('chr',':').split(':')
  if kind!='SNP':continue
  k=(chrom,int(pos),frozenset((r['a1'],r['a2'])));lookup.setdefault(k,[]).append(r)
 except (ValueError,KeyError):continue
base=Path('work/stage4_magma_reference/g1000_eur');found=[];anchors=[]
with open(str(base)+'.bim') as f:
 for i,line in enumerate(f):
  c,rs,_,p,a,b=line.split();p=int(p)
  if c!='15':continue
  if p==74134883 and frozenset((a,b))==frozenset(('C','T')):anchors.append((i,rs,a,b))
  for x in lookup.get((c,p,frozenset((a,b))),[]):found.append((i,rs,p,a,b,x))
assert len(anchors)==1
n=sum(1 for _ in open(str(base)+'.fam'));step=(n+3)//4;bed=np.memmap(str(base)+'.bed',mode='r',dtype=np.uint8);assert bytes(bed[:3])==bytes([0x6c,0x1b,0x01])
lut=np.array([2,np.nan,1,0]);shift=np.array([0,2,4,6],dtype=np.uint8)
def decode(i):return lut[((bed[3+i*step:3+(i+1)*step,None]>>shift)&3).ravel()[:n]]
a=decode(anchors[0][0]);out=[]
for i,rs,p,aa,bb,x in found:
 b=decode(i);ok=np.isfinite(a)&np.isfinite(b)
 if ok.sum()<400 or np.std(a[ok])==0 or np.std(b[ok])==0:continue
 rr=np.corrcoef(a[ok],b[ok])[0,1];out.append(dict(reference_id=rs,godmc_rsid=x['rsid'],chr=15,pos=p,alleles=aa+'/'+bb,reference_pairwise_n=int(ok.sum()),r2=float(rr*rr),godmc_p=x['pval'],godmc_cis=x['cistrans']))
out.sort(key=lambda a:a['r2'],reverse=True)
with open('outputs/stage6_loxl1_godmc_selected_reference_ld.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print('godmc rows',len(rows),'position allele matched',len(found),'evaluated',len(out),'maximum',out[0] if out else None)
