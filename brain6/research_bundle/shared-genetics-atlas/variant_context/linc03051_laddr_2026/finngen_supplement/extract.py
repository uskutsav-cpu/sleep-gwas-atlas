#!/usr/bin/env python3
"""Read publisher's selected FinnGen molecular-colocalization workbook."""
import csv,hashlib,json
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent
src=Path('/Volumes/Brain6Work/laddr-public-supplement-20261028/mmc4.xlsx')
archive=src.parent/'supplements.zip'
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
w=openpyxl.load_workbook(src,read_only=True,data_only=True);s=w['Table S3']
it=s.iter_rows(values_only=True);header=next(it)
assert header==('tissue','gene_id','modality','phenotype_id','finngen_trait','finngen_region','cs_overlap','clpp','clpa','top_shared_variant')
gene=[];traits=set();n=0
for row in it:
 n+=1;traits.add(row[4])
 if row[1]=='ENSG00000239268':gene.append(row)
assert n==63741 and len(gene)==12
with (p/'gene_rows.tsv').open('w',newline='') as f:
 writer=csv.writer(f,delimiter='\t');writer.writerow(header);writer.writerows(gene)
with (p.parent/'sleep_hits.tsv').open() as f:sleep=list(csv.DictReader(f,delimiter='\t'))
matched=[]
for row in gene:
 for x in sleep:
  if row[0]==x['TISSUE'] and row[3]==x['ID'] and row[2]=='latent_residual' and x['model_set']=='gtex-rddp':
   matched.append({'tissue':row[0],'feature':row[3],'finngen_trait':row[4],'finngen_clpp':row[7],
    'finngen_clpa':row[8],'finngen_top_shared_variant':row[9],
    'insomnia_twas_p':x['TWAS.P'],'insomnia_coloc_pp4':x['COLOC.PP4']})
assert len(matched)==2
with (p/'matched_sleep_models.tsv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=matched[0].keys(),delimiter='\t');writer.writeheader();writer.writerows(matched)
receipt={'source_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13361010/supplementaryFiles',
 'archive_path':str(archive),'archive_sha256':sha(archive),'workbook_path':str(src),'workbook_sha256':sha(src),
 'sheet':'Table S3','rows_read':n,'unique_reported_traits':len(traits),'gene_rows':len(gene),
 'gene_endpoint_counts':{x:sum(r[4]==x for r in gene) for x in sorted({r[4] for r in gene})},
 'reported_trait_strings_containing_adhd_or_hyperkinetic':[x for x in sorted(traits) if any(y in str(x).upper() for y in ('ADHD','HYPERKIN'))],
 'sleep_model_matches':len(matched),'gene_rows_sha256':sha(p/'gene_rows.tsv'),'matched_sha256':sha(p/'matched_sleep_models.tsv'),
 'note':'Selected published CLPP >0.01 table; absent entries are not evidence of no effect or of being untested.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
