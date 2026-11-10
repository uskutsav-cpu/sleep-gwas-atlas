import csv,json,pathlib,pysam,numpy as np
rows=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv'))); freq={r['variant']:r for r in csv.DictReader(open('outputs/22_candidate_cohort_allele_frequency_check.csv'))};panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));pops={r['sample']:r['pop'] for r in panel if r['super_pop']=='EUR'};base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};res=[]
for x in rows:
 chrom=x['GRCh37_chr_assumed'];pos=int(x['GRCh37_pos_assumed']);a=x['counted_allele'];b=x['other_allele'];rs=x['variant']
 if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
 f=files[chrom];h=[v for v in f.fetch(chrom,pos-1,pos) if v.id==rs and v.pos==pos and len(v.alleles)==2 and set(v.alleles)=={a,b}];assert len(h)==1
 v=h[0];ix=v.alleles.index(a);d={}
 for s in sorted(set(pops)&set(f.header.samples)):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):d[s]=sum(g==ix for g in gt)
 assert len(d)==503
 u=float(freq[rs]['UKB_frequency_assuming_source_EAF_UKB_is_source_EA']);t=float(freq[rs]['TwinGene_frequency_assuming_source_EAF_TG_is_source_EA']);out={'variant':rs,'counted_allele':a,'ukb_frequency':u,'tg_reported_frequency':t}
 for pop in ('FIN','GBR'):
  xg=[d[s] for s in d if pops[s]==pop];freqp=float(np.mean(xg)/2);out.update({pop+'_n':len(xg),pop+'_counted_allele_frequency':freqp,pop+'_tg_reported_distance':abs(t-freqp),pop+'_tg_complement_distance':abs(1-t-freqp),pop+'_evidence_margin':abs(t-freqp)-abs(1-t-freqp)})
 res.append(out)
out=pathlib.Path('outputs/22_candidate_frequency_reference_subgroups.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=res[0]);w.writeheader();w.writerows(res)
summary={p:{'n':res[0][p+'_n'],'closer_complement':sum(x[p+'_evidence_margin']>0 for x in res),'closer_reported':sum(x[p+'_evidence_margin']<0 for x in res),'absolute_margin_under_0.05':[x['variant'] for x in res if abs(x[p+'_evidence_margin'])<.05]} for p in ('FIN','GBR')}
pathlib.Path('outputs/22_candidate_frequency_reference_subgroups.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
