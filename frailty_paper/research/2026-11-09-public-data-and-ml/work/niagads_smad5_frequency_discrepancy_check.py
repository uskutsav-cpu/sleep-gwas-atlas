import csv, json, hashlib, pathlib, pysam
from collections import Counter
p=pathlib.Path('outputs/stage6_niagads_smad5_exc_Mak_FI_aligned.csv')
rows=[r for r in csv.DictReader(p.open()) if abs(float(r['qtl_alt_af'])-float(r['gwas_af_alt']))>0.1]
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')); eur={x['sample'] for x in panel if x['super_pop']=='EUR'}
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
f=pysam.VariantFile(url); result=[]
for r in rows:
 chrom,pos38,ref38,alt38=r['variant_GRCh38'].split(':'); pos=int(r['gwas_pos37']); found=[]
 for v in f.fetch('5',pos-1,pos):
  if v.pos!=pos or len(v.alleles)!=2 or set(v.alleles)!={ref38,alt38}:continue
  idx=v.alleles.index(alt38); dos=[]
  for s in sorted(eur.intersection(f.header.samples)):
   g=v.samples[s].get('GT')
   if g and len(g)==2 and all(i is not None for i in g):dos.append(sum(i==idx for i in g))
  if dos:found.append({'rsid':v.id,'reference_alleles':list(v.alleles),'n':len(dos),'alt38_frequency':sum(dos)/(2*len(dos))})
 result.append({**r,'reference_matches':found})
out=pathlib.Path('outputs/stage6_niagads_smad5_frequency_discrepancy.json');out.write_text(json.dumps({'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(pathlib.Path('work/1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'rows':result,'note':'Descriptive EUR reference frequency and allele-match check only; does not verify strand, cohort ancestry, source coding or causality. The remote VCF full file was not checksum verified.'},indent=2)+'\n')
print(json.dumps(result,indent=2))
