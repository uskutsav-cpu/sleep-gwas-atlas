import csv,json,hashlib,openpyxl
from pathlib import Path
files={'GEPITEL.zip':'5a01fa78c023ffb96c6eaaab4c6aad32','GWAS_LMI.xlsx':'c58c6d9b4427fc9b16c223fe4bc1faae'}
for name,md in files.items():
 h=hashlib.md5()
 with open('work/gepitel/'+name,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 assert h.hexdigest()==md,(name,h.hexdigest())
a=list(csv.DictReader(open('outputs/integrated_15_candidate_review.csv')));b=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')))
ids={r['variant'] for r in a+b};positions={}
for label,rows in [('15',a),('22',b)]:
 for r in rows:
  for build in ['37','38']:
   chr=r['GRCh37_chr'] if label=='15' else (r['GRCh37_chr_assumed'] if build=='37' else r['GRCh38_chr']); pos=r['GRCh'+build+'_pos'] if label=='15' else r['GRCh'+build+'_pos_assumed'] if build=='37' else r['GRCh38_pos'];positions.setdefault((str(chr),str(pos)),[]).append((label,build,r['variant']))
found={'bim_id':[],'bim_coordinate':[],'lmi_id':[],'lmi_coordinate':[]}
for line in open('work/gepitel/GEPITEL.bim'):
 r=line.split()
 if r[1] in ids:found['bim_id'].append(r[1])
 if (r[0],r[3]) in positions:found['bim_coordinate'].extend(positions[(r[0],r[3])])
w=openpyxl.load_workbook('work/gepitel/GWAS_LMI.xlsx',read_only=True,data_only=True);s=w['dataset'];n=0
for r in s.iter_rows(min_row=2,values_only=True):
 n+=1;chrom,pos=str(r[0]),str(r[2]);rs=r[1]
 if rs in ids:found['lmi_id'].append({'variant':rs,'chr':chrom,'position':pos,'p':str(r[6])})
 if (chrom,pos) in positions:found['lmi_coordinate'].extend(positions[(chrom,pos)])
result={'file_md5_verified':files,'fam_n':192,'fam_phenotypes_present':0,'bim_variants':654027,'lmi_rows':n,'findings':found,'limitations':'Unspecified genome build, unknown genotype imputation coverage and no linked phenotype; no absence of association can be inferred.'}
Path('outputs/stage6_gepitel_public_data_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
