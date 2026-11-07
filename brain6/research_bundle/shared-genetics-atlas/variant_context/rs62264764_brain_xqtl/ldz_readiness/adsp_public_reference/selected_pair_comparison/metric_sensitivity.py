#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
import numpy as np,pysam
p=Path(__file__).resolve().parent;atlas=p.parents[4];prior=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot';panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
def read(q):
 with q.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
receipt=json.loads((prior/'ld_receipt.json').read_text());assert sha(panel)==receipt['panel_sha256'] and sha(prior/'variants.tsv')==receipt['variants_sha256']
ids=sorted(x['sample'] for x in read(panel) if x['super_pop']=='EUR');assert len(ids)==503
pairs=read(p/'bootstrap.tsv');assert len(pairs)==3
v={r['variant']:r for r in read(prior/'variants.tsv')};keys=sorted({x[k] for x in pairs for k in ('variant1','variant2')})
vcf=pysam.VariantFile(receipt['reference_url']); assert set(ids)<=set(vcf.header.samples);hap={}
for key in keys:
 x=v[key];pos=int(x['pos37']);records=[q for q in vcf.fetch('3',pos-1,pos) if q.pos==pos and q.id and x['rsid'] in q.id.split(';') and x['ref'] in q.alleles and x['alt'] in q.alleles];assert len(records)==1
 q=records[0];ri=q.alleles.index(x['ref']);ai=q.alleles.index(x['alt']);h=[]
 for ident in ids:
  call=q.samples[ident];gt=call.get('GT');assert gt and len(gt)==2 and call.phased and all(g in (ri,ai) for g in gt);h.extend(int(g==ai) for g in gt)
 hap[key]=np.asarray(h,dtype=float).reshape(503,2)
out=[]
for x in pairs:
 a,b=hap[x['variant1']],hap[x['variant2']];dos=float(np.corrcoef(a.sum(1),b.sum(1))[0,1]**2);h=float(np.corrcoef(a.ravel(),b.ravel())[0,1]**2);assert abs(dos-float(x['reference_r2']))<1e-9
 out.append({'variant1':x['variant1'],'variant2':x['variant2'],'dosage_r2':dos,'haplotype_r2':h,'absolute_metric_difference':abs(dos-h),'adsp_r2':x['adsp_r2']})
with (p/'metric_sensitivity.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(p/'metric_receipt.json').write_text(json.dumps({'panel_sha256':sha(panel),'input_bootstrap_sha256':sha(p/'bootstrap.tsv'),'individuals':503,'phased_haplotypes':1006,'reference_url':receipt['reference_url'],'whole_vcf_and_index_hashed':False,'selection':'postselected top three between-reference discrepancies'},indent=2)+'\n');print((p/'metric_sensitivity.tsv').read_text())
