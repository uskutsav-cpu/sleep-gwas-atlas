import pysam,csv,numpy as np,json
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
with open('work/1kg_phase3_panel.tsv') as f:eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
c=list(csv.DictReader(open('outputs/candidate_risk_gene_review_15.csv')))
pairs=[('rs113851554','rs78206187'),('rs35881094','rs17550821'),('rs11152363','rs1942262'),('rs11152363','rs660010'),('rs1942262','rs660010')]
lookup={r['variant']:(r['chr_GRCh37'],int(r['position_GRCh37']),r['counted_allele'],r['other_allele']) for r in c};lookup['rs660010']=('18',53310969,'G','C')
needed={v for pair in pairs for v in pair}; files={}; data={};status={};meta={}
for rs in sorted(needed):
 chrom,pos,a,b=lookup[rs]
 try:
  if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
  v=files[chrom];samples=sorted(eur & set(v.header.samples));matches=[x for x in v.fetch(chrom,pos-1,pos) if x.id==rs and x.pos==pos]
  if len(matches)!=1:raise ValueError('matching VCF records: '+str(len(matches)))
  x=matches[0];meta[rs]={'alleles':list(x.alleles),'target':a,'other':b,'chr':chrom,'pos':pos,'samples':len(samples)}
  if set(x.alleles)!={a,b}:raise ValueError('reference VCF allele pair not matching study contrast '+str(x.alleles))
  ix=x.alleles.index(a); d={}
  for sample in samples:
   g=x.samples[sample].get('GT')
   if g and len(g)==2 and all(i is not None for i in g) and x.samples[sample].phased:d[sample]=[int(q==ix) for q in g]
  data[rs]=d;status[rs]='matched_exact_VCF_allele_pair';print(rs,status[rs],meta[rs]['alleles'],flush=True)
 except Exception as e:status[rs]='unavailable: '+str(e);print(rs,status[rs],flush=True)
res=[]
for a,b in pairs:
 r={'variant_1':a,'counted_allele_1':lookup[a][2],'variant_2':b,'counted_allele_2':lookup[b][2],'status':'','n_individuals':None,'n_haplotypes':None,'r':None,'r2':None}
 if a not in data or b not in data:r['status']='unavailable: '+status[a]+'; '+status[b]
 else:
  s=sorted(data[a].keys() & data[b].keys());x=np.array([x for i in s for x in data[a][i]]);y=np.array([x for i in s for x in data[b][i]])
  if x.std()==0 or y.std()==0:r['status']='zero allele variance'
  else:
   cor=float(np.corrcoef(x,y)[0,1]);r.update(status='matched biallelic VCF, phased EUR',n_individuals=len(s),n_haplotypes=len(x),r=cor,r2=cor*cor)
 res.append(r)
 print(r,flush=True)
with open('outputs/candidate_regional_allele_specific_ld.json','w') as f:json.dump({'reference':'1000 Genomes phase3 v5a, UCSC hg19','population':'EUR','variant_status':status,'variant_metadata':meta,'pairs':res,'limitation':'Reference VCF contrast only; study multiallelic coding still needs verification; no clinical or causal inference'},f,indent=2)
