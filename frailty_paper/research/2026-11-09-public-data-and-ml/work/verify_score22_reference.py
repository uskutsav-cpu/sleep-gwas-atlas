import csv,json,hashlib,importlib.util,pathlib,pysam,numpy as np
src=pathlib.Path('outputs/score_22_research.py');spec=importlib.util.spec_from_file_location('score22',src);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
weights=list(csv.DictReader(open('outputs/unvalidated_22_variant_score_prototype.csv')));panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'}
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={}; data={}
for row in weights:
 chrom=row['chr'];pos=int(row['position_GRCh37_assumed']);a=row['counted_allele'];b=row['other_allele'];rs=row['variant']
 if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
 f=files[chrom];hits=[v for v in f.fetch(chrom,pos-1,pos) if v.id==rs and v.pos==pos and len(v.alleles)==2 and set(v.alleles)=={a,b}];assert len(hits)==1
 v=hits[0];idx=v.alleles.index(a);data[rs]={}
 for s in eur&set(f.header.samples):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):data[rs][s]=sum(g==idx for g in gt)
samples=sorted(set.intersection(*(set(x) for x in data.values())));assert len(samples)==503
obs=[]
for s in samples:
 for i,r in enumerate(weights):
  d=data[r['variant']][s];reverse=(i%2==0)
  obs.append({'sample_id':s,'variant':r['variant'],'dosage_allele':r['other_allele'] if reverse else r['counted_allele'],'other_allele':r['counted_allele'] if reverse else r['other_allele'],'dosage':str(2-d if reverse else d)})
score=module.calculate(weights,obs,allow_verified_reversal=True); assert len(score)==503
arr=np.array([[v['FI_raw_sum'],v['FFS_raw_sum'],v['UKB_insomnia_raw_sum']] for v in score]);corr=np.corrcoef(arr,rowvar=False);old=np.asarray(json.load(open('outputs/reference_22_vs_15_score_diagnostic.json'))['new22_cross_trait_score_correlations']);delta=float(np.max(np.abs(corr-old)));assert delta<1e-12
result={'samples':503,'n_variants':len(weights),'reversed_input_contrasts':11,'correlation_matrix':corr.tolist(),'maximum_absolute_difference_from_prior':delta,'script_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'weights_sha256':hashlib.sha256(pathlib.Path('outputs/unvalidated_22_variant_score_prototype.csv').read_bytes()).hexdigest(),'warning':'Checks arithmetic and allele reversal against the prior genotype-only reference calculation; no outcome data and no prediction validation. No sample genotypes or scores persisted.'};pathlib.Path('outputs/score22_reference_implementation_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
