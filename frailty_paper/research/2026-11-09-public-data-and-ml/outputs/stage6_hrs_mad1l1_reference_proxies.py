#!/usr/bin/env python3
"""Exhaustive local screen against the 32 published HRS leads, not full GWAS."""
import csv,json
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent; out=root/'outputs'
def read(name):
    with (out/name).open(newline='') as f:return list(csv.DictReader(f))
leads=read('stage6_mekli2018_published_fi_leads.csv');score22=read('unvalidated_22_shared_direction_score.csv');score15=read('unvalidated_15_variant_score_prototype.csv')
pairs=[]
for label,score,ch,pos in [('22',score22,'chr_GRCh37_assumed','position_GRCh37_assumed'),('15',score15,'chr','position_GRCh37_inferred')]:
    for s in score:
        for h in leads:
            if s[ch]==h['chr'] and abs(int(s[pos])-int(h['bp']))<=1000000:
                pairs.append((label,s,h))
assert len(leads)==32 and len(pairs)==6 and {p[1]['variant'] for p in pairs}=={'rs28528997'}
with (root/'work/1kg_phase3_panel.tsv').open() as f:names=sorted(r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR')
assert len(names)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url)
def dosages(pos,rs,count=None):
    matches=[r for r in vcf.fetch('7',pos-1,pos) if r.pos==pos and r.id==rs]
    if len(matches)!=1:return None,{'status':'rsid_position_not_unique','count':len(matches)}
    r=matches[0]
    if len(r.alleles)!=2:return None,{'status':'not_biallelic','alleles':list(r.alleles)}
    count=count or r.alleles[1]
    if count not in r.alleles:return None,{'status':'counted_allele_mismatch','alleles':list(r.alleles)}
    idx=r.alleles.index(count);v=[]
    for sample in names:
        gt=r.samples[sample].get('GT');v.append(sum(a==idx for a in gt) if gt is not None and len(gt)==2 and all(a is not None for a in gt) else np.nan)
    return np.array(v,dtype=float),{'status':'matched','alleles':list(r.alleles),'counted':count}
anchor,anchor_info=dosages(2074122,'rs28528997','A');assert anchor_info['status']=='matched'
rows=[]
for label,s,h in pairs:
    v,info=dosages(int(h['bp']),h['snp']);row={'candidate':s['variant'],'lead':h['snp'],'distance_bp':int(h['bp'])-int(s['position_GRCh37_assumed']),'lead_hrs_a1':h['a1_hrs'],'hrs_beta':h['beta_hrs'],'hrs_p':h['p_hrs'],'elsa_beta':h['beta_elsa'],'elsa_p':h['p_elsa'],'reference_status':info['status'],'reference_alleles':'/'.join(info.get('alleles',[])),'n_complete':None,'r_squared':None}
    if v is not None:
        keep=np.isfinite(anchor)&np.isfinite(v)
        if keep.sum()>2 and np.std(v[keep])>0:
            row['n_complete']=int(keep.sum());row['r_squared']=float(np.corrcoef(anchor[keep],v[keep])[0,1]**2)
    rows.append(row)
with (out/'stage6_hrs_mad1l1_reference_proxies.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'stage6_hrs_mad1l1_reference_proxies.json').write_text(json.dumps({'input_scope':'32 selected HRS rows; 15 and 22 prototype candidates; all pairs <=1Mb','source':url,'reference_population':'1000 Genomes EUR phase 3','anchor':anchor_info,'pairs':rows,'limitations':'Reference r2 is not HRS LD; lead allele strand and exact pair not confirmed in published HRS table; selected leads not genome-wide coverage; no score validation or individual outcomes.'},indent=2)+'\n')
print(json.dumps(rows,indent=2))
