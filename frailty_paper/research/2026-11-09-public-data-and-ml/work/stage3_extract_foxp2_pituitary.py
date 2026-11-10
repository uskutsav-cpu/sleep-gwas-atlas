import pysam,pathlib,csv,hashlib,json,time
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000301/QTD000301.all.tsv.gz';idx=pathlib.Path('work/QTD000301.all.tsv.gz.tbi'); start,end=113086327,115086327;gene='ENSG00000128573';cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split();f=pysam.TabixFile(url,index=str(idx));rows=0;keep=[];started=time.time()
for line in f.fetch('7',start-1,end):
 rows+=1;v=line.split('\t')
 if len(v)!=len(cols):raise ValueError('unexpected columns')
 if v[16].split('.')[0]==gene:keep.append(dict(zip(cols,v)))
 if rows%1000000==0:print('read',rows,'kept',len(keep),flush=True)
p=pathlib.Path('work/stage3_QTD000301_FOXP2_2mb.tsv')
with p.open('w',newline='') as o:w=csv.DictWriter(o,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(keep)
receipt={'url':url,'index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'catalogue_metadata_release':7,'dataset':'QTD000301','phenotype':'pituitary gene expression, n=237 per release metadata','gene':gene,'build':'GRCh38','region_inclusive':[start,end],'region_rows_read':rows,'gene_rows':len(keep),'elapsed_seconds':time.time()-started,'gene_table_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'caveat':'Remote indexed subset, full source and index checksum not externally verified; limited window; not a completed colocalization'}
pathlib.Path('outputs/stage3_catalogue_foxp2_pituitary_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
