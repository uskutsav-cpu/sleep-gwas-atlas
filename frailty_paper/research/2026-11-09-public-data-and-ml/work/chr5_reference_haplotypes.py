import csv,pysam,numpy as np,json,hashlib,math
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
panel='work/1kg_phase3_panel.tsv'
with open(panel) as f:
 eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
v=pysam.VariantFile(url); samples=sorted(eur&set(v.header.samples));print('EUR samples',len(samples),flush=True)
variants={'rs1592757':(103889998,'C','G'),'rs7703746':(103907337,'G','A'),'rs323509':(104082179,'A','C')}; vals={}; meta={}
for rs,(pos,target,other) in variants.items():
 records=[r for r in v.fetch('5',pos-1,pos) if r.id==rs]
 if len(records)!=1:raise ValueError((rs,len(records)))
 r=records[0];a=r.alleles
 if set(a)!={target,other}:raise ValueError((rs,a))
 targetidx=a.index(target);out=[]
 for s in samples:
  call=r.samples[s];gt=call.get('GT');
  if not gt or len(gt)!=2 or any(i is None for i in gt) or not call.phased:out.extend([np.nan,np.nan])
  else:out.extend([int(i==targetidx) for i in gt])
 vals[rs]=np.array(out,dtype=float);meta[rs]={'pos':pos,'alleles':list(a),'counted':target,'other':other,'counted_freq':float(np.nanmean(vals[rs]))}
 print(rs,meta[rs],flush=True)
res=[]
for a,b in [('rs323509','rs7703746'),('rs323509','rs1592757'),('rs7703746','rs1592757')]:
 x,y=vals[a],vals[b];m=np.isfinite(x)&np.isfinite(y);x,y=x[m],y[m];r=float(np.corrcoef(x,y)[0,1]);
 counts={'00':int(np.sum((x==0)&(y==0))),'01':int(np.sum((x==0)&(y==1))),'10':int(np.sum((x==1)&(y==0))),'11':int(np.sum((x==1)&(y==1)))}
 res.append({'variant_1':a,'counted_allele_1':meta[a]['counted'],'variant_2':b,'counted_allele_2':meta[b]['counted'],'sample_count':len(samples),'complete_haplotypes':len(x),'r':r,'r_squared':r*r,'counts':counts})
print(json.dumps(res,indent=2))
with open('outputs/chr5_allele_specific_reference_ld.json','w') as f:json.dump({'vcf_url':url,'panel_url':'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel','panel_sha256':hashlib.sha256(open(panel,'rb').read()).hexdigest(),'sample_population':'EUR','variants':meta,'pairs':res,'warning':'Reference haplotypes only; no study-specific phase or causal colocalization established'},f,indent=2)
