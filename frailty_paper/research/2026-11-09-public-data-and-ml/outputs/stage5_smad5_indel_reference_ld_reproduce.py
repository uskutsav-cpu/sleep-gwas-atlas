import csv,collections,json,pathlib,hashlib,datetime,math
import numpy as np,pysam
p=pathlib.Path('outputs');url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';index='work/stage5_1000g_chr5_phase3.tbi';t=pysam.TabixFile(url,index=index);hdr=[h for h in t.header if h.startswith('#CHROM')];assert len(hdr)==1;samples=hdr[0].split('\t')[9:]
panel={x['sample']:x for x in csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')};groups={name:[i for i,s in enumerate(samples) if s in panel and pred(panel[s])] for name,pred in [('EUR',lambda x:x['super_pop']=='EUR'),('FIN',lambda x:x['pop']=='FIN')]};assert len(groups['EUR'])==503 and len(groups['FIN'])==99
cs=[x for x in csv.DictReader(open(p/'stage5_credible_qtl_standardized_deduplicated.csv')) if x['dataset_id']=='QTD000285' and x['gene_id']=='ENSG00000113658' and len(x['effect_allele'])!=len(x['other_allele'])];maps={x['variant_GRCh38']:x for x in csv.DictReader(open(p/'stage5_indel_hg19_mapping.csv'))};gen={};inventory=[]
for x in cs:
 k=x['variant_GRCh38'];m=maps[k];pos=int(m['hg19_position']);rows=[]
 for l in t.fetch('5',pos-1,pos):
  z=l.split('\t');alts=z[4].split(',')
  if int(z[1])==pos and z[3]==m['hg19_ref'] and m['hg19_alt'] in alts:rows.append(z)
 if len(rows)!=1:inventory.append({'variant_GRCh38':k,'status':'no unique exact REF/ALT record','ref_rsid':'','reference_alt_af_EUR':'','reference_alt_af_FIN':''});continue
 z=rows[0];al=z[4].split(',').index(m['hg19_alt'])+1;fmt=z[8].split(':');gi=fmt.index('GT');dos=[]
 for sample in z[9:]:
  gt=sample.split(':')[gi].replace('|','/').split('/');dos.append(sum(int(a)==al for a in gt) if len(gt)==2 and all(a.isdigit() for a in gt) else float('nan'))
 v=np.array(dos,dtype=float);gen[k]=v
 def af(idx):a=v[idx];return float(np.nanmean(a)/2) if np.isfinite(a).any() else ''
 inventory.append({'variant_GRCh38':k,'status':'unique exact REF/ALT record','ref_rsid':z[2],'reference_alt_af_EUR':af(groups['EUR']),'reference_alt_af_FIN':af(groups['FIN'])})
with open(p/'stage5_smad5_indel_reference_eur_dosages.csv','w',newline='') as df:
 w=csv.writer(df);w.writerow(['sample','population','variant_GRCh38','alt_dosage'])
 for k,v in gen.items():
  for i in groups['EUR']:
   w.writerow([samples[i],panel[samples[i]]['pop'],k,'' if not math.isfinite(v[i]) else int(v[i])])
keys=list(gen);pairs=[]
for group,idx in groups.items():
 for i in range(len(keys)):
  for j in range(i+1,len(keys)):
   a=gen[keys[i]][idx];b=gen[keys[j]][idx];valid=np.isfinite(a)&np.isfinite(b);av=a[valid];bv=b[valid];corr=float(np.corrcoef(av,bv)[0,1]) if len(av)>2 and np.std(av)>0 and np.std(bv)>0 else float('nan');pairs.append({'population':group,'variant_1':keys[i],'variant_2':keys[j],'n_nonmissing':int(valid.sum()),'r_alt_dosage':corr if math.isfinite(corr) else '', 'r_squared':corr*corr if math.isfinite(corr) else '', 'status':'estimated' if math.isfinite(corr) else 'monomorphic or insufficient data'})
for fn,rs in [('stage5_smad5_indel_reference_inventory.csv',inventory),('stage5_smad5_indel_reference_pairwise_ld.csv',pairs)]:
 with open(p/fn,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_vcf':url,'source_index_sha256':hashlib.sha256(pathlib.Path(index).read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(pathlib.Path('work/1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'sample_counts':{k:len(v) for k,v in groups.items()},'credible_indels':len(cs),'exact_vcf_alleles':len(gen),'missing':len(cs)-len(gen),'pairwise_rows':len(pairs),'method':'unphased ALT dosage Pearson r for pairwise complete public reference genotypes','limitations':'small and ancestry-specific reference; 1000G EUR includes FIN; missing allele not proven absent, and 7-indel LD is not a complete locus LD matrix'};(p/'stage5_smad5_indel_reference_ld_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
