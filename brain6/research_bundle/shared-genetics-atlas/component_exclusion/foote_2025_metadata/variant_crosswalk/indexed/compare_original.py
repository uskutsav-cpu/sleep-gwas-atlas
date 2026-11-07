import csv,pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent
with (p/'GCST90624046_anchors.tsv').open() as f:h={x['rsid_requested']:x for x in csv.DictReader(f,delimiter='\t')}
out=[]
for rs in ('rs62264764','rs817509','rs2431108','rs7703746'):
 with (p/f'{rs}_original_row.tsv').open() as f:o=next(csv.DictReader(f,delimiter='\t'))
 x=h[rs];assert x['status']=='UNIQUE_MATCH' and o['rs_id']==rs and {o['effect_allele'],o['other_allele']}=={x['effect_allele'],x['other_allele']}
 flip=1 if o['effect_allele']==x['effect_allele'] else -1
 aligned=float(o['beta'])*flip;af=float(o['effect_allele_frequency']);af=af if flip==1 else 1-af
 assert abs(aligned-float(x['beta']))<1e-12 and abs(af-float(x['effect_allele_frequency']))<1e-12
 for k in ('standard_error','p_value','Q_metric','Q_metric_pvalue'):assert abs(float(o[k])-float(x[k]))<1e-12*max(1,abs(float(o[k])))
 out.append({'rsid':rs,'original_EA':o['effect_allele'],'harmonized_EA':x['effect_allele'],'beta_alignment':flip,'aligned_beta':aligned,'harmonized_beta':x['beta'],'aligned_eaf':af,'harmonized_eaf':x['effect_allele_frequency'],'p':x['p_value'],'Q_p':x['Q_metric_pvalue'],'comparison':'PASS'})
with (p/'original_harmonized_comparison.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
print('Four matching variant rows: allele-adjusted beta, frequency, SE, p and Q fields agree.')
