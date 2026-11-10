import csv,json,hashlib,pathlib,pysam,numpy as np
orig=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv'))); weights={r['variant']:r for r in csv.DictReader(open('outputs/22_fi_cohort_weight_sensitivity.csv'))};panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')); eur={r['sample'] for r in panel if r['super_pop']=='EUR'};pops={r['sample']:r['pop'] for r in panel};base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};G=[]
for x in orig:
 c=x['GRCh37_chr_assumed'];p=int(x['GRCh37_pos_assumed']);a=x['counted_allele'];b=x['other_allele'];assert weights[x['variant']]['counted_allele']==a
 if c not in files:files[c]=pysam.VariantFile(base.format(chrom=c))
 f=files[c];hits=[v for v in f.fetch(c,p-1,p) if v.id==x['variant'] and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={a,b}];assert len(hits)==1
 v=hits[0];ix=v.alleles.index(a);d={}
 for s in sorted(eur & set(f.header.samples)):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):d[s]=sum(g==ix for g in gt)
 G.append(d)
samples=sorted(set.intersection(*(set(d) for d in G)));assert len(samples)==503
X=np.array([[d[s] for d in G] for s in samples],float); U=np.array([float(weights[x['variant']]['ukb_beta']) for x in orig]); T=np.array([float(weights[x['variant']]['twingene_beta_inferred']) for x in orig]);scores=np.column_stack((X@U,X@T));overall=float(np.corrcoef(scores,rowvar=False)[0,1]);assert abs(overall-(-0.4547968355928501))<1e-12
res=[]
for i,x in enumerate(orig):
 sub=scores-np.column_stack((X[:,i]*U[i],X[:,i]*T[i]));r=float(np.corrcoef(sub,rowvar=False)[0,1]);res.append({'excluded_variant':x['variant'],'ukb_twingene_score_correlation':r,'difference_from_full':r-overall})
out=pathlib.Path('outputs/22_fi_score_cohort_leave_one_out.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=res[0]);w.writeheader();w.writerows(res)
subs={}
for p in sorted(set(pops[s] for s in samples)):
 ix=[i for i,s in enumerate(samples) if pops[s]==p]; subs[p]={'n':len(ix),'correlation':float(np.corrcoef(scores[ix,:],rowvar=False)[0,1])}
exclude={'rs12536335','rs28528997','rs17312861','rs11636381'};ix=[i for i,x in enumerate(orig) if x['variant'] not in exclude];drop4=float(np.corrcoef(np.column_stack((X[:,ix]@U[ix],X[:,ix]@T[ix])),rowvar=False)[0,1]);assert len(ix)==18
j={'full_correlation':overall,'n':len(samples),'leave_one_out_min':min(res,key=lambda x:x['ukb_twingene_score_correlation']),'leave_one_out_max':max(res,key=lambda x:x['ukb_twingene_score_correlation']),'drop_four_selected_heterogeneity_candidates_correlation':drop4,'drop_four_variants':sorted(exclude),'subgroup_descriptive':subs,'output_csv_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'limitations':'Post hoc; no phenotype data, not predictive validity, inferred TwinGene sign adapter, small subgroups, selection and ancestry.'};pathlib.Path('outputs/22_fi_score_cohort_leave_one_out.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
