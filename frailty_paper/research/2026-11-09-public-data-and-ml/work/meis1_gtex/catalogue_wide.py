import hashlib,json
from pathlib import Path
import pysam
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000166/QTD000166.all.tsv.gz'
out=Path('work/meis1_gtex/QTD000166_MEIS1AS3_65400000_67500000.tsv')
f=pysam.TabixFile(url);total=0;kept=0
with out.open('w') as o:
 for line in f.fetch('2',65399999,67500000):
  total+=1
  if line.split('\t',1)[0]=='ENSG00000226819':o.write(line+'\n');kept+=1
receipt={'url':url,'region_GRCh38':'2:65400000-67500000','gene_id':'ENSG00000226819','total_region_rows':total,'matching_gene_rows':kept,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'metadata_QTD000166_release7_sample_size':209,'caveat':'wide region based on NCBI coordinates; exact gene annotation and full cis coverage remain unverified; not independent of GTEx v10 donors'}
Path('outputs/stage6_meis1_catalogue_cerebellum_wide_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
