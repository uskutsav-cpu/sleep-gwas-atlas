#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import numpy as np,pysam
p=Path(__file__).resolve().parent;atlas=p.parents[4];prior=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot';panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
receipt=json.loads((prior/'ld_receipt.json').read_text());assert sha(panel)==receipt['panel_sha256'] and sha(prior/'variants.tsv')==receipt['variants_sha256']
ids=sorted(x['sample'] for x in read(panel) if x['super_pop']=='EUR');assert len(ids)==503
pairs=sorted([r for r in read(p/'pair_comparison.tsv') if r['status']=='matched'],key=lambda x:-float(x['absolute_difference']))[:3]
v={r['variant']:r for r in read(prior/'variants.tsv')};keys=sorted({x[k] for x in pairs for k in ('variant1','variant2')})
url=receipt['reference_url'];vcf=pysam.VariantFile(url);assert set(ids)<=set(vcf.header.samples)
G={}
for key in keys:
 x=v[key];pos=int(x['pos37']);found=[q for q in vcf.fetch('3',pos-1,pos) if q.pos==pos and q.id and x['rsid'] in q.id.split(';') and x['ref'] in q.alleles and x['alt'] in q.alleles];assert len(found)==1
 q=found[0];ri=q.alleles.index(x['ref']);ai=q.alleles.index(x['alt']);arr=[]
 for ident in ids:
  gt=q.samples[ident].get('GT');assert gt and len(gt)==2 and all(g in (ri,ai) for g in gt);arr.append(sum(g==ai for g in gt))
 G[key]=np.asarray(arr,dtype=float)
rng=np.random.default_rng(20261028);resamples=rng.integers(0,len(ids),size=(5000,len(ids)));out=[]
for x in pairs:
 a,b=G[x['variant1']],G[x['variant2']];r=float(np.corrcoef(a,b)[0,1]**2);assert abs(r-float(x['kgp_r2']))<1e-9
 aa=a[resamples];bb=b[resamples];aa-=aa.mean(axis=1,keepdims=True);bb-=bb.mean(axis=1,keepdims=True);den=np.sqrt((aa*aa).sum(axis=1)*(bb*bb).sum(axis=1));good=den>0
 rr=((aa[good]*bb[good]).sum(axis=1)/den[good])**2;lo,hi=np.quantile(rr,[.025,.975]);adsp=float(x['adsp_r2'])
 out.append({'variant1':x['variant1'],'variant2':x['variant2'],'reference_r2':r,'adsp_r2':adsp,'bootstrap_low_2_5':lo,'bootstrap_high_97_5':hi,'valid_replicates':len(rr),'monomorphic_replicates':5000-len(rr),'fraction_at_or_below_adsp':float((rr<=adsp).mean())})
with (p/'bootstrap.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(p/'bootstrap_receipt.json').write_text(json.dumps({'vcf_url':url,'panel_sha256':sha(panel),'pair_comparison_sha256':sha(p/'pair_comparison.tsv'),'selection':'three largest previously observed discrepancies','n':503,'replicates':5000,'seed':20261028,'remote_vcf_not_whole_file_hashed':True,'limitation':'Within-reference resampling only; postselected; no sample or LD-statistic equivalence.'},indent=2)+'\n');print((p/'bootstrap.tsv').read_text())
