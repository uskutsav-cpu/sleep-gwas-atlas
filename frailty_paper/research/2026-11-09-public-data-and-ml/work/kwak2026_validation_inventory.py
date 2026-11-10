import csv,openpyxl,hashlib,pathlib,json
p=pathlib.Path('outputs');source=pathlib.Path('work/kwak2026_supplement.xlsx');w=openpyxl.load_workbook(source,read_only=True,data_only=True)
sets={'prototype_15':{x['variant'] for x in csv.DictReader((p/'unvalidated_15_variant_score_prototype.csv').open())},'expanded_22':{x['variant'] for x in csv.DictReader((p/'expanded_22_candidate_risk_gene_table.csv').open())},'song_124':{x['variant_id'] for x in csv.DictReader((p/'published_insomnia_fi_124_cross_trait_candidates.csv').open())}}
rows=[]
for sheet,trait in [('SData2','FRAIL scale'),('SData3','CHS index'),('SData4','Frailty index')]:
 s=w[sheet]
 for i,row in enumerate(s.iter_rows(min_row=4,values_only=True),4):
  if not row[0] or not str(row[0]).strip().startswith('rs'):continue
  x={'sheet':sheet,'source_row':i,'phenotype':trait,'rsid':str(row[0]).strip(),'chromosome':row[1],'bp_reported_build_unverified':str(row[2]).strip(),'source_gwas_catalog_allele':row[3],'source_gwas_catalog_beta':row[5],'source_gwas_catalog_p':row[6],'KoGES_allele_reported':row[7],'KoGES_beta_reported':row[9],'KoGES_se_reported':row[10],'KoGES_p_reported':row[11]}
  if len(row)>12:x.update(KNHANES_allele_reported=row[12],KNHANES_beta_reported=row[14],KNHANES_se_reported=row[15],KNHANES_p_reported=row[16],Meta_allele_reported=row[17],Meta_beta_reported=row[19],Meta_se_reported=row[20],Meta_p_reported=row[21])
  x.update({k+'_exact_rsid_match':x['rsid'] in v for k,v in sets.items()});rows.append(x)
q=p/'kwak2026_supplementary_candidate_inventory.csv';keys=list(dict.fromkeys(k for r in rows for k in r))
with q.open('w',newline='') as f:wri=csv.DictWriter(f,fieldnames=keys);wri.writeheader();wri.writerows(rows)
print('sheet_rows', {k:sum(x['sheet']==k for x in rows) for k in ('SData2','SData3','SData4')},'overlap',{k:sorted({r['rsid'] for r in rows if r[k+'_exact_rsid_match']}) for k in sets})
print('source_sha256',hashlib.sha256(source.read_bytes()).hexdigest(),'output_sha256',hashlib.sha256(q.read_bytes()).hexdigest())
print('docx_sha256',hashlib.sha256(pathlib.Path('work/kwak2026_supplement.docx').read_bytes()).hexdigest())
