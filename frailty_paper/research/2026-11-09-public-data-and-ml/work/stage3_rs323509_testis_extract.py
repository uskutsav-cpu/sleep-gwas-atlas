import pysam,pathlib,csv,hashlib,json,time,collections
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000340/QTD000340.cc.tsv.gz';idx=pathlib.Path('work/QTD000340.cc.tsv.gz.tbi');start,end=104246478,105246478;gene='ENSG00000251574';cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split();f=pysam.TabixFile(url,index=str(idx));rows=0;keep=[];started=time.time()
for line in f.fetch('5',start-1,end):
 rows+=1;v=line.split('\t');assert len(v)==len(cols)
 if v[16].split('.')[0]==gene:keep.append(dict(zip(cols,v)))
p=pathlib.Path('work/stage3_QTD000340_rs323509_testis.tsv')
with p.open('w',newline='') as o:w=csv.DictWriter(o,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(keep)
receipt={'url':url,'index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'catalogue_release':7,'dataset':'QTD000340','phenotype':'testis leafcutter splicing, 322 samples per release metadata','gene':gene,'build':'GRCh38','region_inclusive':[start,end],'region_rows_read':rows,'gene_rows':len(keep),'event_counts':dict(collections.Counter(x['molecular_trait_id'] for x in keep)),'elapsed_seconds':time.time()-started,'table_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'caveat':'filtered cc file selected molecular traits only; remote index and full source not independently authenticated; v8 Catalogue cluster labels may differ from v10 GTEx portal'}
pathlib.Path('outputs/stage3_rs323509_testis_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
