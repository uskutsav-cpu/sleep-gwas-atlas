import pathlib,csv,json,hashlib,collections,time,pysam,math
p=pathlib.Path(__file__).resolve().parent
base='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/'
queries=[('Brain_Nucleus_accumbens_basal_ganglia','ENSG00000176349'),('Brain_Caudate_basal_ganglia','ENSG00000122687')]
columns=['variant','r2','pvalue','molecular_trait_object_id','molecular_trait_id','maf','gene_id','median_tpm','beta','se','an','ac','chromosome','position','ref','alt','type','rsid']
summary=[]
for tissue,gene in queries:
 url=base+tissue+'.tsv.gz';reader=pysam.TabixFile(url);dest=p/(tissue+'_'+gene+'.tsv');c=collections.Counter();snps=set();duplicates=0
 with dest.open('w',newline='') as f:
  w=csv.DictWriter(f,columns,delimiter='\t');w.writeheader()
  for line in reader.fetch('7',1534486,2534487):
   c['region_rows_all_genes']+=1;row=line.split('\t');assert len(row)==len(columns);x=dict(zip(columns,row))
   if x['molecular_trait_id']!=gene:continue
   assert x['gene_id']==gene and x['chromosome']=='7' and 1534487<=int(x['position'])<=2534487
   assert x['variant']==f"chr7_{x['position']}_{x['ref']}_{x['alt']}"
   try:
    pv=float(x['pvalue']);beta=float(x['beta']);se=float(x['se'])
    if not (math.isfinite(pv) and 0<=pv<=1 and all(map(math.isfinite,(beta,se))) and se>0):raise ValueError('out-of-range or invalid')
    c['valid_statistic_rows']+=1
   except ValueError:c['invalid_statistic_rows']+=1
   c['gene_rows']+=1
   if x['variant'] in snps:duplicates+=1
   snps.add(x['variant']);w.writerow(x)
 summary.append({'tissue':tissue,'gene':gene,'url':url,'region_grch38':'7:1534487-2534487','rows_scanned':c['region_rows_all_genes'],'gene_rows':c['gene_rows'],'valid_statistic_rows':c['valid_statistic_rows'],'invalid_statistic_rows':c['invalid_statistic_rows'],'unique_variant_keys':len(snps),'duplicate_variant_rows':duplicates,'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'file':dest.name});print(summary[-1],flush=True);time.sleep(2)
(p/'receipt.json').write_text(json.dumps({'queries':summary,'source_note':'remote indexed GTEx v8 imported into eQTL Catalogue; whole data and index not hashed','format':'https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md','retrieved_date':'2026-10-27'},indent=2)+'\n')
