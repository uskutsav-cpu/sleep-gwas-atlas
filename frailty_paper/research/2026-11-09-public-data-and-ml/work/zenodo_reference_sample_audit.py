import csv,json,pathlib,hashlib
import pysam
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'};assert len(eur)==503
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};out=[]
for r in csv.DictReader(open('work/zenodo_genomewide_systematic_rows.csv')):
 c=r['chromosome'];pos=int(r['base_pair_location']);e=r['effect_allele'];o=r['other_allele'];record={'row_number':r['row_number'],'variant':r['rsid'],'chr':c,'pos':pos,'source_eaf':r['effect_allele_frequency'],'status':'','paired_calls':'','reference_eaf':'','absolute_difference':'','within_4dp':'','error':''}
 try:
  if c not in files:files[c]=pysam.VariantFile(base.format(chrom=c))
  f=files[c];match=[x for x in f.fetch(c,pos-1,pos) if x.id==r['rsid'] and len(x.alleles)==2 and set(x.alleles)=={e,o}]
  if len(match)!=1:record['status']='no unique exact biallelic reference match'
  else:
   v=match[0];idx=v.alleles.index(e);calls=[]
   for name in eur:
    gt=v.samples[name].get('GT')
    if gt and len(gt)==2 and all(x is not None for x in gt):calls.append(sum(x==idx for x in gt))
   if len(calls)!=503:record['status']='incomplete EUR calls';record['paired_calls']=len(calls)
   else:
    q=sum(calls)/(2*len(calls));diff=abs(float(r['effect_allele_frequency'])-q);record.update(status='matched',paired_calls=len(calls),reference_eaf=q,absolute_difference=diff,within_4dp=diff<=.0000500000001)
 except Exception as ex:record['status']='reference retrieval error';record['error']=str(ex)[:180]
 out.append(record);print(r['row_number'],record['status'],flush=True)
for f in files.values():f.close()
with open('outputs/stage6_zenodo_systematic_reference_audit.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
valid=[r for r in out if r['status']=='matched']; summary={'method':'systematic rows at 123456 + k*300000 in file order, selected before looking at association or frequencies; not random','selected_rows':len(out),'unique_exact_reference_matches':len(valid),'matching_within_0_00005':sum(r['within_4dp'] for r in valid),'reference_retrieval_errors':sum(r['status']=='reference retrieval error' for r in out),'max_difference':max([r['absolute_difference'] for r in valid],default=None),'selected_input_sha256':hashlib.sha256(pathlib.Path('work/zenodo_genomewide_systematic_rows.csv').read_bytes()).hexdigest(),'caveat':'selected file rows with unique matching biallelic variants and 503 complete EUR calls; not proof of source or error in effects'}
pathlib.Path('outputs/stage6_zenodo_systematic_reference_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
