import csv,gzip,json,pathlib,collections
w=[x for x in csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')) if x['candidate']];receipts=json.load(open('outputs/stage3_remaining_catalogue_extract_receipt.json'));res={x['candidate']:x for x in receipts['items']};bydataset=collections.defaultdict(list)
for x in w:bydataset[x['dataset_r7']].append(x)
selected=[];summ=[]
for d,items in bydataset.items():
 path=pathlib.Path('work')/(d+'.credible_sets.tsv.gz');geneids={x['gene_id'] for x in items};records=[]
 if path.exists():
  with gzip.open(path,'rt') as f:
   for x in csv.DictReader(f,delimiter='\t'):
    if x['gene_id'].split('.')[0] in geneids:records.append(x)
 for item in items:
  c=item['candidate'];subset=[x for x in records if x['gene_id'].split('.')[0]==item['gene_id']];groups=collections.defaultdict(dict)
  for x in subset:
   k=x['variant'];key=x['cs_id'];pip=float(x['pip']);assert k not in groups[key] or abs(groups[key][k]-pip)<1e-12;groups[key][k]=pip;selected.append({'candidate':c,'dataset':d,'gene_symbol':item['gene_symbol'],**x})
  ext=res[c];lo,hi=ext['region'];unique={v for z in groups.values() for v in z};within={v for v in unique if lo<=int(v.split('_')[1])<=hi};cvar='chr'+item['marker_grch38'].split(':')[0].removeprefix('chr')+'_'+item['marker_grch38'].split(':')[1]+'_'+item['candidate_ref']+'_'+item['candidate_alt'];
  summ.append({'candidate':c,'dataset':d,'gene_id':item['gene_id'],'gene_symbol':item['gene_symbol'],'tissue':item['tissue_v10'],'region_gene_rows':ext.get('gene_rows',0),'exact_candidate_rows':ext.get('candidate_rows',0),'credible_set_count':len(groups),'credible_set_sizes':';'.join(k+':'+str(len(v)) for k,v in sorted(groups.items())),'credible_set_unique_variants':len(unique),'credible_set_variants_in_candidate_window':len(within),'candidate_in_any_set':cvar in unique,'source_cs_file_exists':path.exists(),'status':('no gene rows in regional QTL source' if ext.get('gene_rows',0)==0 else 'regional QTL available; molecular set reported' if groups else 'regional QTL available; no molecular credible set reported'),'limitation':'molecular model only; no disease causality; v8-v10 related data; QTL candidate window may truncate fine-mapping window'})
with open('outputs/stage3_remaining_molecular_credible_sets.csv','w',newline='') as f:
 if selected:wri=csv.DictWriter(f,fieldnames=selected[0]);wri.writeheader();wri.writerows(selected)
with open('outputs/stage3_remaining_molecular_finemapping_summary.csv','w',newline='') as f:wri=csv.DictWriter(f,fieldnames=summ[0]);wri.writeheader();wri.writerows(summ)
for x in summ:print(x['candidate'],x['gene_symbol'],'qtl rows',x['region_gene_rows'],'marker rows',x['exact_candidate_rows'],'CS',x['credible_set_count'],'variants',x['credible_set_unique_variants'],'candidate',x['candidate_in_any_set'])
