import csv,math,pysam,numpy as np,json,hashlib,pathlib
out=pathlib.Path('outputs');new=list(csv.DictReader((out/'expanded_22_candidate_risk_gene_table.csv').open()));old=list(csv.DictReader((out/'unvalidated_15_variant_score_prototype.csv').open()));assert len(new)==22 and len(old)==15
fields=['variant','chr','position_GRCh37_assumed','counted_allele','other_allele','FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight','nearest_gene_body_positional','reference_alleles','status']
rows=[]
for x in new:rows.append(dict(variant=x['variant'],chr=x['GRCh37_chr_assumed'],position_GRCh37_assumed=x['GRCh37_pos_assumed'],counted_allele=x['counted_allele'],other_allele=x['other_allele'],FI_beta_weight=x['FI_beta_aligned'],FFS_beta_weight_inferred=x['FFS_beta_aligned'],UKB_insomnia_beta_weight=x['UKB_insomnia_beta_aligned'],nearest_gene_body_positional=x['nearest_gene_body_positional'],reference_alleles=x['reference_alleles'],status='Postselected unvalidated trait-specific illustration; no individual risk prediction; verify alleles, LD and weights'))
proto=out/'unvalidated_22_variant_score_prototype.csv'
with proto.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
newmap={r['variant']:r for r in rows};oldmap={r['variant']:r for r in old};union=sorted(set(newmap)|set(oldmap));assert len(union)==33
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'};pops={r['sample']:r['pop'] for r in panel};base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};data={};checks=[]
for v in union:
 r=newmap.get(v) or oldmap[v]; chrom=r['chr'];pos=int(r.get('position_GRCh37_assumed') or r['position_GRCh37_inferred']);a=r['counted_allele'];b=r['other_allele']
 if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
 f=files[chrom];hits=[x for x in f.fetch(chrom,pos-1,pos) if x.id==v and x.pos==pos and len(x.alleles)==2 and set(x.alleles)=={a,b}]
 assert len(hits)==1,(v,len(hits));x=hits[0];ix=x.alleles.index(a);data[v]={}
 for sample in sorted(eur&set(f.header.samples)):
  gt=x.samples[sample].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):data[v][sample]=sum(g==ix for g in gt)
 checks.append({'variant':v,'counted_allele':a,'other_allele':b,'calls':len(data[v])});print(v,'matched',len(data[v]),flush=True)
samples=sorted(set.intersection(*(set(x) for x in data.values())));assert len(samples)==503

def calc(rows, cols):
 arr=np.zeros((len(samples),len(cols)))
 for r in rows:
  v=r['variant']; ref=newmap.get(v) or oldmap[v]; raw=np.array([data[v][s] for s in samples],float)
  if r['counted_allele']!=ref['counted_allele']:
   assert r['counted_allele']==ref['other_allele'] and r['other_allele']==ref['counted_allele'];raw=2-raw
  for j,key in enumerate(cols):arr[:,j]+=raw*float(r[key])
 return arr
A=calc(old,['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']);B=calc(rows,['FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight'])
prev=json.load(open(out/'reference_15_score_diagnostic.json'));assert np.allclose(np.corrcoef(A,rowvar=False),prev['score_pearson_correlations'],rtol=0,atol=1e-12)
combo=np.corrcoef(np.column_stack((A,B)),rowvar=False); vals={}
for i,name in enumerate(['FI','FFS','UKB_insomnia']): vals[name]=float(combo[i,i+3])
subs={}
for pop in sorted(set(pops[s] for s in samples)):
 ix=[i for i,s in enumerate(samples) if pops[s]==pop];subs[pop]={'n':len(ix),'same_trait_correlations':{name:float(np.corrcoef(A[ix,j],B[ix,j])[0,1]) for j,name in enumerate(['FI','FFS','UKB_insomnia'])}}
res={'status':'genotype-only technical diagnostic, no outcomes or prospective risk','n_reference_eur_samples':len(samples),'n_old_variants':len(old),'n_new_variants':len(rows),'n_overlap':len(set(newmap)&set(oldmap)),'n_union':len(union),'same_trait_correlations_old15_vs_new22':vals,'new22_cross_trait_score_correlations':np.corrcoef(B,rowvar=False).tolist(),'subgroups':subs,'variant_checks':checks,'old_score_file_sha256':hashlib.sha256((out/'unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'new_score_file_sha256':hashlib.sha256(proto.read_bytes()).hexdigest(),'warning':'Scores post-selected, GWAS cohorts overlap or overlap unknown, possible multiallelic encoding, FFS orientation inferred, no phenotype outcomes, not causal or clinical.'}
(out/'reference_22_vs_15_score_diagnostic.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='variant_checks'},indent=2))
