import pysam,gzip,hashlib,json
from pathlib import Path
acc='GCST90824107';u=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/{acc}/harmonised/{acc}.h.tsv.gz';idx=Path('work')/(acc+'.h.tsv.gz.tbi');f=pysam.TabixFile(u,index=str(idx));out=Path('work/stage6_portas_nrem_foxp2_region.tsv.gz');n=0
with gzip.open(out,'wt') as t:
 for row in f.fetch('7',113903000-1,114903000):
  t.write(row+'\n');n+=1
print('rows',n,flush=True);Path('outputs/stage6_portas_nrem_foxp2_region_receipt.json').write_text(json.dumps({'url':u,'index_md5':hashlib.md5(idx.read_bytes()).hexdigest(),'chromosome':'7','start_GRCh38':113903000,'end_GRCh38':114903000,'rows':n,'retrieved_file_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'full_source_checksum_verified':False,'note':'Tabix extracted subset of harmonised GRCh38 records, no header, only sex-combined NREM'},indent=2)+'\n')
