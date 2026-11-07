import pathlib,csv,json,hashlib,collections,time,pysam,math
p=pathlib.Path(__file__).resolve().parent;base='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/'
tissues=('Testis','Brain_Frontal_Cortex_BA9');genes=('ENSG00000251574',)
cols=['variant','r2','pvalue','molecular_trait_object_id','molecular_trait_id','maf','gene_id','median_tpm','beta','se','an','ac','chromosome','position','ref','alt','type','rsid'];summary=[]
for tissue in tissues:
 url=base+tissue+'.tsv.gz';reader=pysam.TabixFile(url);handles={};writers={};counts={g:collections.Counter() for g in genes};seen={g:set() for g in genes};scanned=0
 for g in genes:
  f=(p/(tissue+'_'+g+'.tsv')).open('w',newline='');handles[g]=f;writers[g]=csv.DictWriter(f,cols,delimiter='\t');writers[g].writeheader()
 try:
  for line in reader.fetch('5',104112266,105112267):
   scanned+=1;values=line.split('\t');assert len(values)==len(cols);x=dict(zip(cols,values));g=x['molecular_trait_id']
   if g not in genes:continue
   assert x['gene_id']==g and x['chromosome']=='5' and 104112267<=int(x['position'])<=105112267 and x['variant']==f"chr5_{x['position']}_{x['ref']}_{x['alt']}"
   counts[g]['rows']+=1
   if x['variant'] in seen[g]:counts[g]['duplicate_variant_rows']+=1
   seen[g].add(x['variant'])
   try:
    pv=float(x['pvalue']);b=float(x['beta']);se=float(x['se']);maf=float(x['maf'])
    if not (math.isfinite(pv) and 0<=pv<=1 and math.isfinite(b) and math.isfinite(se) and se>0 and math.isfinite(maf) and 0<maf<=.5):raise ValueError()
    counts[g]['valid_stats']+=1
   except ValueError:counts[g]['invalid_stats']+=1
   writers[g].writerow(x)
 finally:
  for f in handles.values():f.close()
 for g in genes:
  f=p/(tissue+'_'+g+'.tsv');summary.append({'tissue':tissue,'gene':g,'url':url,'region_grch38':'5:104112267-105112267','rows_scanned_all_genes':scanned,'counts':dict(counts[g]),'unique_variant_keys':len(seen[g]),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'file':f.name})
 print(tissue,[(g,dict(counts[g])) for g in genes],flush=True);time.sleep(2)
(p/'qtl_receipt.json').write_text(json.dumps({'queries':summary,'source_note':'indexed remote eQTL Catalogue imported GTEx v8; full source and index not hashed','retrieval_date':'2026-10-27'},indent=2)+'\n')
