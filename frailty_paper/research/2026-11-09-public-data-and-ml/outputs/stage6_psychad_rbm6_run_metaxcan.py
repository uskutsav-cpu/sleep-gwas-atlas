#!/usr/bin/env python3
from pathlib import Path
import csv,subprocess,json
root=Path(__file__).resolve().parents[1];out=root/'outputs';work=root/'work';traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_insomnia','FinnGen_HFRS_R12'];rows=[];receipts={}
for name in traits:
 result=work/f'stage6_psychad_rbm6_{name}_result.csv';cmd=['python3',str(work/'MetaXcan/software/SPrediXcan.py'),'--model_db_path',str(work/'stage6_psychad_rbm6_model.db'),'--covariance',str(work/'stage6_psychad_rbm6_cov.txt.gz'),'--gwas_file',str(work/f'stage6_rbm6_{name}_gwas.tsv'),'--snp_column','rsid','--effect_allele_column','effect_allele','--non_effect_allele_column','other_allele','--beta_column','beta','--se_column','se','--output_file',str(result),'--throw','--overwrite'];r=subprocess.run(cmd,cwd=root,text=True,capture_output=True);receipts[name]={'returncode':r.returncode,'stderr':r.stderr[-6000:]};r.check_returncode();rr=list(csv.DictReader(result.open()));assert len(rr)==1 and rr[0]['gene']=='ENSG00000004534';rows.append({'trait':name,**rr[0]})
with open(out/'stage6_psychad_rbm6_twas_exploratory.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'stage6_psychad_rbm6_metaxcan_run_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
for x in rows:print(x['trait'],'Z',x['zscore'],'P',x['pvalue'],'SNPs',x['n_snps_used'])
