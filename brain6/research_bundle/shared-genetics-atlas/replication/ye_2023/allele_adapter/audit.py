import csv,hashlib,json
from pathlib import Path
import openpyxl
here=Path(__file__).resolve().parent
source=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
wb=openpyxl.load_workbook(here/'source.xlsx',read_only=True,data_only=True)
rows=list(wb[' ST1. COJO independent loci'].iter_rows(values_only=True)); hdr=rows[0]; leads=[dict(zip(hdr,row)) for row in rows[1:] if row[0] is not None]; ids={str(x['SNP']) for x in leads}; found={i:[] for i in ids}; n=0
with source.open() as f:
 reader=csv.DictReader(f,delimiter='\t')
 for r in reader:
  n+=1
  if r['SNP'] in ids: found[r['SNP']].append(r)
results=[]
for lead in leads:
 match=[r for r in found[str(lead['SNP'])] if r['CHROM']==str(lead['Chr']) and r['POS']==str(lead['bp'])]
 rec={'snp':lead['SNP'],'chr':lead['Chr'],'bp':lead['bp'],'refA':lead['refA'],'b':lead['b'],'se':lead['se'],'matches':len(match)}
 if len(match)==1:
  r=match[0]; rec.update({'A1':r['A1'],'A2':r['A2'],'BETA':r['BETA'],'SE':r['SE'],'allele_equal':r['A1']==lead['refA'],'beta_equal':float(r['BETA'])==float(lead['b']),'se_equal':float(r['SE'])==float(lead['se'])})
 results.append(rec)
with (here/'comparisons.tsv').open('w',newline='') as f:
 keys=['snp','chr','bp','refA','b','se','matches','A1','A2','BETA','SE','allele_equal','beta_equal','se_equal']; w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(results)
receipt={'workbook_sha256':sha(here/'source.xlsx'),'previous_workbook_sha256':'cc95a22e0c8731473d29c72951ee3b887bc9344f87bbb01f73dc4ee8709282e9','source_path':str(source),'source_sha256':sha(source),'source_rows':n,'lead_rows':len(leads),'unique_matches':sum(x['matches']==1 for x in results),'allele_equal':sum(x.get('allele_equal',False) for x in results),'beta_equal':sum(x.get('beta_equal',False) for x in results),'se_equal':sum(x.get('se_equal',False) for x in results)}
(here/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
