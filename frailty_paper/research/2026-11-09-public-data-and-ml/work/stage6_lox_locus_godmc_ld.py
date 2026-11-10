import json,csv,hashlib
import numpy as np
from pathlib import Path
p=Path('work/stage4_magma_reference/g1000_eur');ids={'rs11636381','rs4077284'}; records={}
with open(str(p)+'.bim') as f:
 for i,l in enumerate(f):
  c,rs,_,bp,a,b=l.split()
  if rs in ids:records.setdefault(rs,[]).append((i,c,int(bp),a,b))
assert all(len(records.get(x,[]))==1 for x in ids)
ni=sum(1 for _ in open(str(p)+'.fam'));step=(ni+3)//4;bed=np.memmap(str(p)+'.bed',mode='r',dtype='uint8');assert bytes(bed[:3])==bytes([0x6c,0x1b,0x01])
def decode(k):
 x=records[k][0];codes=((bed[3+x[0]*step:3+(x[0]+1)*step,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).ravel()[:ni];return np.array([2,np.nan,1,0])[codes]
a=decode('rs11636381');b=decode('rs4077284');ok=np.isfinite(a)&np.isfinite(b);rr=np.corrcoef(a[ok],b[ok])[0,1]
j=json.load(open('work/godmc_cpg_cg12594244.json'));rows=j['assoc_meta'];sel=[r for r in rows if r.get('rsid')=='rs11636381'];lead=[r for r in rows if r.get('rsid')=='rs4077284'];assert len(lead)==1
r={'reference':'MAGMA 1000 Genomes EUR; GRCh37','candidate':records['rs11636381'][0][1:],'godmc_mqtl':records['rs4077284'][0][1:],'pairwise_n':int(ok.sum()),'r2':float(rr*rr),'signed_r_bim_a1':float(rr),'godmc_cpg_rows':len(rows),'target_exact_rows':len(sel),'godmc_lead':lead[0],'transport_note':'HTTP source: no TLS authentication; reference LD is not target study LD. GoDMC table is preselected, not full nonsignificant results.'}
Path('outputs/stage6_loxl1_godmc_ld.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
