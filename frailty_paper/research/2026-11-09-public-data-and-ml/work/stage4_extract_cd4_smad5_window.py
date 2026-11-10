import pysam,pathlib,csv,hashlib,json,time
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000612/QTD000612.all.tsv.gz';idx=pathlib.Path('work/QTD000612.all.tsv.gz.tbi'); start,end=135702653,136702653;gene='ENSG00000113658';cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split();f=pysam.TabixFile(url,index=str(idx));rows=0;keep=[];started=time.time()
for line in f.fetch('5',start-1,end):
 rows+=1;v=line.split('\t')
 if len(v)!=len(cols):raise ValueError('unexpected columns')
 if v[16].split('.')[0]==gene:keep.append(dict(zip(cols,v)))
 if rows%1000000==0:print('read',rows,'kept',len(keep),flush=True)
p=pathlib.Path('work/stage4_QTD000612_SMAD5_1mb.tsv')
with p.open('w',newline='') as o:w=csv.DictWriter(o,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(keep)
receipt={'url':url,'index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'catalogue_metadata_release':7,'dataset':'QTD000612','phenotype':'naive CD4 T cell pseudobulk gene expression, n=949 per release metadata','gene':gene,'build':'GRCh38','region_inclusive':[start,end],'region_rows_read':rows,'gene_rows':len(keep),'elapsed_seconds':time.time()-started,'gene_table_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'caveat':'Remote indexed subset, full source and index checksum not externally verified; limited window; not a completed colocalization'}
pathlib.Path('outputs/stage4_catalogue_cd4_smad5_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
