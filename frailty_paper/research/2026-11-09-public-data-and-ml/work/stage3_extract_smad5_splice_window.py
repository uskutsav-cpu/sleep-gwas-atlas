import pysam,pathlib,csv,hashlib,json,time
url='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000285/QTD000285.cc.tsv.gz';idx=pathlib.Path('work/QTD000285.cc.tsv.gz.tbi'); start,end=135702653,136702653;event='5:136132962:136153592:clu_31306_+';cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split();f=pysam.TabixFile(url,index=str(idx));rows=0;keep=[];started=time.time()
for line in f.fetch('5',start-1,end):
 rows+=1;v=line.split('\t');assert len(v)==len(cols)
 if v[0]==event:keep.append(dict(zip(cols,v)))
p=pathlib.Path('work/stage3_QTD000285_SMAD5_splice_1mb.tsv')
with p.open('w',newline='') as o:w=csv.DictWriter(o,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(keep)
receipt={'url':url,'index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'catalogue_metadata_release':7,'dataset':'QTD000285','phenotype':'skeletal muscle leafcutter splicing, 702 samples per release metadata','gene':'ENSG00000113658','event':event,'build':'GRCh38','region_inclusive':[start,end],'region_rows_read':rows,'event_rows':len(keep),'elapsed_seconds':time.time()-started,'event_table_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'caveat':'Filtered cc file retains selected strongest molecular traits from fine-mapped signals; exact event not comprehensive across all SMAD5 splice events. Remote indexed subset, full source and index checksum not independently verified.'}
pathlib.Path('outputs/stage3_catalogue_smad5_splicing_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
