import hashlib,json
from pathlib import Path
import pysam
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000166/QTD000166.all.tsv.gz'
out=Path('work/meis1_gtex/QTD000166_MEIS1AS3_66400000_66530000.tsv')
f=pysam.TabixFile(url);total=0;kept=0
with out.open('w') as o:
 for line in f.fetch('2',66399999,66530000):
  total+=1
  if line.split('\t',1)[0]=='ENSG00000226819':o.write(line+'\n');kept+=1
receipt={'url':url,'region_GRCh38':'2:66400000-66530000','gene_id':'ENSG00000226819','total_region_rows':total,'matching_gene_rows':kept,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'metadata_QTD000166_release7_sample_size':209,'caveat':'region subset from live tabix, not full cis window; gene transcripts and source version require audit; not independent of GTEx v10 donors'}
Path('outputs/stage6_meis1_catalogue_cerebellum_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
