from pathlib import Path
import pysam,csv,json,hashlib,time
O=Path('outputs');W=Path('work');u='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz';idx=W/'QTD000171.all.tsv.gz.tbi';f=pysam.TabixFile(u,index=str(idx));cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split();rows=[];n=0;t=time.time()
for line in f.fetch('3',49714152-1,50711834):
 x=line.split('\t');assert len(x)==len(cols);n+=1
 if x[16].split('.')[0]=='ENSG00000004534':rows.append(dict(zip(cols,x)))
 if n%500000==0:print('region rows',n,'gene rows',len(rows),flush=True)
p=O/'stage6_rbm6_gtex_cortex_region.tsv';
with p.open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
rec={'source_url':u,'index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'index_checksum_verified_against_source':False,'whole_source_checksum_verified':False,'gene':'ENSG00000004534','region_GRCh38':'chr3:49714152-50711834','n_regional_rows':n,'n_gene_rows':len(rows),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'dataset':'QTD000171 GTEx v8 bulk brain cortex, catalog metadata N=205','elapsed_seconds':time.time()-t};(O/'stage6_rbm6_gtex_cortex_extract_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2),flush=True)
