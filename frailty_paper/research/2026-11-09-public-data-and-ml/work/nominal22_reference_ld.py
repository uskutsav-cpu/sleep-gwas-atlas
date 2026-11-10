import csv,json,pysam,numpy as np,itertools,pathlib,hashlib
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
rows=[r for r in csv.DictReader(open('outputs/song124_mvp_hfrs_joint_screen.csv')) if r.get('joint_p_max') and float(r['joint_p_max'])<.05];assert len(rows)==22
ids={r['sample'] for r in csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t') if r['super_pop']=='EUR'};files={};data={};status={}
for r in rows:
 v=r['variant'];c=r['song_chr'];p=int(r['song_position_reported']);a=r['HFRS_effect_allele'];b=r['HFRS_other_allele']
 try:
  if c not in files:files[c]=pysam.VariantFile(base.format(chrom=c))
  f=files[c];samples=sorted(ids&set(f.header.samples)); hits=[x for x in f.fetch(c,p-1,p) if x.id==v and x.pos==p and len(x.alleles)==2 and set(x.alleles)=={a,b}]
  if len(hits)!=1:raise ValueError(f'exact ID position biallelic pair matches: {len(hits)}')
  x=hits[0]; ix=x.alleles.index(a);calls={}
  for sample in samples:
   g=x.samples[sample].get('GT')
   if g and len(g)==2 and all(k is not None for k in g) and x.samples[sample].phased:calls[sample]=[int(k==ix) for k in g]
  data[v]=calls;status[v]={'status':'matched','n_called':len(calls),'ref_alleles':list(x.alleles)}
 except Exception as e:status[v]={'status':'unavailable: '+str(e)}
 print(v,status[v]['status'],flush=True)
pairs=[]
for x,y in itertools.combinations(rows,2):
 if x['song_chr']!=y['song_chr'] or abs(int(x['song_position_reported'])-int(y['song_position_reported']))>1000000:continue
 a,b=x['variant'],y['variant'];r={'variant_1':a,'variant_2':b,'chr':x['song_chr'],'distance_bp':abs(int(x['song_position_reported'])-int(y['song_position_reported'])),'counted_allele_1':x['HFRS_effect_allele'],'counted_allele_2':y['HFRS_effect_allele'],'status':''}
 if a not in data or b not in data:r['status']='unavailable: '+status[a]['status']+'; '+status[b]['status']
 else:
  common=sorted(data[a].keys() & data[b].keys());aa=np.array([h for s in common for h in data[a][s]]);bb=np.array([h for s in common for h in data[b][s]])
  if len(aa)==0 or aa.std()==0 or bb.std()==0:r['status']='zero variance or missing'
  else:
   cor=float(np.corrcoef(aa,bb)[0,1]);r.update(status='matched phased EUR',n_individuals=len(common),n_haplotypes=len(aa),r=cor,r2=cor*cor)
 pairs.append(r)
out=pathlib.Path('outputs/song124_22_nominal_reference_ld.json');out.write_text(json.dumps({'reference':'1000G phase3 v5a indexed UCSC GRCh37','population':'EUR, pooled 503 reference individuals; no outcome data','source_candidates_sha256':hashlib.sha256(pathlib.Path('outputs/song124_mvp_hfrs_joint_screen.csv').read_bytes()).hexdigest(),'variant_status':status,'pairs_same_chromosome_at_most_1Mb':pairs,'warning':'Conditional on exact biallelic VCF representation; r2 does not prove conditional independence, shared causal effects, or score validity'},indent=2)+'\n')
print('pairs',json.dumps(pairs),flush=True)
