import csv,gzip,pathlib,collections,hashlib,json
sources=[('SMIM32 expression',pathlib.Path('work/QTD000281.credible_sets.tsv.gz'),'ENSG00000271824',None),('SMAD5 selected splice event',pathlib.Path('work/QTD000285.credible_sets.tsv.gz'),'ENSG00000113658','5:136132962:136153592:clu_31306_+')];summary=[];out=[];receipts=[]
for tag,p,gene,event in sources:
 groups=collections.defaultdict(list);total=0
 with gzip.open(p,'rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   total+=1
   if x['gene_id']==gene and (event is None or x['molecular_trait_id']==event):groups[x['cs_id']].append(x)
 for cs,rows in groups.items():
  for x in rows:out.append({'qtl_phenotype':tag,**x})
  unique={}
  for x in rows:
   if x['variant'] in unique: assert x['pip']==unique[x['variant']]['pip']
   else: unique[x['variant']]=x
  vals=list(unique.values()); assert len(vals)==int(rows[0]['cs_size'])
  candidate=[x for x in vals if x['variant']=='chr5_136202653_A_C'];pip=sum(float(x['pip']) for x in vals)
  outside=sum(not (135702653<=int(x['variant'].split('_')[1])<=136702653) for x in vals)
  summary.append({'qtl_phenotype':tag,'molecular_trait_id':rows[0]['molecular_trait_id'],'cs_id':cs,'reported_cs_size':rows[0]['cs_size'],'observed_rows':len(rows),'unique_variants':len(vals),'duplicates_removed':len(rows)-len(vals),'unique_variant_pip_sum':pip,'cs_min_r2':rows[0]['cs_min_r2'],'region':rows[0]['region'],'variants_outside_previous_1mb_window':outside,'rs7734045_in_cs':bool(candidate),'candidate_pip':candidate[0]['pip'] if candidate else '', 'top_variant':max(vals,key=lambda x:float(x['pip']))['variant'],'top_pip':max(float(x['pip']) for x in vals)})
 receipts.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_rows':total,'selected_cs_count':len(groups),'selected_rows':sum(map(len,groups.values()))})
fields=['qtl_phenotype','molecular_trait_id','gene_id','cs_id','variant','rsid','cs_size','pip','pvalue','beta','se','z','cs_min_r2','region']
for name,records in [('outputs/stage3_skeletal_qtl_credible_set_variants.csv',out),('outputs/stage3_skeletal_qtl_credible_set_summary.csv',summary)]:
 if records:
  with open(name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=(fields if records is out else list(records[0])));w.writeheader();w.writerows(records)
pathlib.Path('outputs/stage3_skeletal_qtl_credible_set_receipt.json').write_text(json.dumps({'sources':receipts,'notes':'Provider credible sets describe molecular phenotypes only; not a GWAS credible set or a GWAS-QTL colocalization; complete full-file hash against provider not verified.'},indent=2)+'\n');print('summary',summary);print('receipts',receipts)
