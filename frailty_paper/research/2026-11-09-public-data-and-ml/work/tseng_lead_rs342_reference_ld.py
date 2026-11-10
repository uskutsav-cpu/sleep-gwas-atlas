import csv,pysam,numpy as np,json,pathlib
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
variant=[('rs34290943',49161660,'T','C'),('rs3774745',50204745,'C','T')]; ids={x['sample'] for x in csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t') if x['super_pop']=='EUR'};f=pysam.VariantFile(url);d={};status=[]
for rs,p,a,b in variant:
 hits=[x for x in f.fetch('3',p-1,p) if x.id==rs and x.pos==p and len(x.alleles)==2 and set(x.alleles)=={a,b}]
 if len(hits)!=1:status.append({'id':rs,'n_exact':len(hits)});continue
 x=hits[0];ix=x.alleles.index(a);calls={}
 for s in ids&set(f.header.samples):
  gt=x.samples[s].get('GT')
  if gt and len(gt)==2 and all(y is not None for y in gt) and x.samples[s].phased:calls[s]=[int(y==ix) for y in gt]
 d[rs]=calls;status.append({'id':rs,'n_exact':1,'n_calls':len(calls),'counted_allele':a})
r=None
if len(d)==2:
 samples=sorted(d[variant[0][0]].keys() & d[variant[1][0]].keys());x=np.array([g for s in samples for g in d[variant[0][0]][s]]);y=np.array([g for s in samples for g in d[variant[1][0]][s]]);assert x.std()>0 and y.std()>0;r=float(np.corrcoef(x,y)[0,1])
out={'reference':url,'population':'1000 Genomes phase 3 EUR-labeled, phased','distance_bp':1043085,'status':status,'r':r,'r2':r*r if r is not None else None,'limitations':'Assumed study build and exact allele matches; reference r2 is not conditional independence, causality or verification of study effect alleles.'};pathlib.Path('outputs/tseng_shared22_nearby_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
