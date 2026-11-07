#!/usr/bin/env python3
import collections,csv,hashlib,json
from pathlib import Path
import numpy as np,pysam
p=Path(__file__).resolve().parent;atlas=p.parents[4];prior=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot';panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
def read(q):
 with q.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
receipt=json.loads((prior/'ld_receipt.json').read_text());assert sha(panel)==receipt['panel_sha256'] and sha(prior/'variants.tsv')==receipt['variants_sha256']
panelrows=sorted((x for x in read(panel) if x['super_pop']=='EUR'),key=lambda x:x['sample']);ids=[x['sample'] for x in panelrows];assert len(ids)==503
pop=collections.defaultdict(list)
for i,x in enumerate(panelrows):pop[x['pop']].append(i)
assert set(pop)=={'CEU','FIN','GBR','IBS','TSI'}
pairs=read(p/'bootstrap.tsv');assert len(pairs)==3;v={r['variant']:r for r in read(prior/'variants.tsv')};keys=sorted({x[k] for x in pairs for k in ('variant1','variant2')})
vcf=pysam.VariantFile(receipt['reference_url']);assert set(ids)<=set(vcf.header.samples);G={}
for key in keys:
 x=v[key];pos=int(x['pos37']);q=[z for z in vcf.fetch('3',pos-1,pos) if z.pos==pos and z.id and x['rsid'] in z.id.split(';') and x['ref'] in z.alleles and x['alt'] in z.alleles];assert len(q)==1
 z=q[0];ri=z.alleles.index(x['ref']);ai=z.alleles.index(x['alt']);a=[]
 for ident in ids:
  gt=z.samples[ident].get('GT');assert gt and len(gt)==2 and all(g in (ri,ai) for g in gt);a.append(sum(g==ai for g in gt))
 G[key]=np.asarray(a,dtype=float)
rng=np.random.default_rng(20261028);out=[]
for x in pairs:
 A,B=G[x['variant1']],G[x['variant2']];assert abs(float(np.corrcoef(A,B)[0,1]**2)-float(x['reference_r2']))<1e-9
 for name,ii in sorted(pop.items()):
  a,b=A[ii],B[ii];r2=float(np.corrcoef(a,b)[0,1]**2) if a.std()>0 and b.std()>0 else float('nan')
  ix=rng.integers(0,len(ii),size=(3000,len(ii)));aa=a[ix];bb=b[ix];aa-=aa.mean(axis=1,keepdims=True);bb-=bb.mean(axis=1,keepdims=True);den=np.sqrt((aa*aa).sum(axis=1)*(bb*bb).sum(axis=1));ok=den>0; vals=((aa[ok]*bb[ok]).sum(axis=1)/den[ok])**2
  lo,hi=np.quantile(vals,[.025,.975]) if len(vals) else (float('nan'),float('nan'))
  out.append({'variant1':x['variant1'],'variant2':x['variant2'],'population':name,'n':len(ii),'alt_frequency_1':float(a.mean()/2),'alt_frequency_2':float(b.mean()/2),'dosage_r2':r2,'adsp_r2':x['adsp_r2'],'bootstrap_low':lo,'bootstrap_high':hi,'valid_bootstrap':len(vals),'monomorphic_bootstrap':3000-len(vals)})
with (p/'population.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(p/'population_receipt.json').write_text(json.dumps({'panel_sha256':sha(panel),'bootstrap_input_sha256':sha(p/'bootstrap.tsv'),'n':len(ids),'population_counts':{k:len(v) for k,v in pop.items()},'n_pairs':len(pairs),'replicates':3000,'seed':20261028,'whole_vcf_and_index_hashed':False,'caution':'Postselected, groups are small; no cross-panel hypothesis test.'},indent=2)+'\n');print((p/'population.tsv').read_text())
