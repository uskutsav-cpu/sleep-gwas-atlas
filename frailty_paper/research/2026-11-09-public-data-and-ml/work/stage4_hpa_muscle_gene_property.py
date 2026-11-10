import csv,hashlib,json,math,pathlib,subprocess
root=pathlib.Path('work/stage4_magma');matrix=root/'matrices/HPA_human_skeletal_muscle_GSE233882.txt';traits=['jansen_insomnia_brown','atkins_fi_brown','mak_fi_brown','ye_ffs_brown'];out=[];receipt=[]
for trait in traits:
 prefix=root/(trait+'_HPA_human_skeletal_muscle_GSE233882');o=pathlib.Path(str(prefix)+'.gsa.out');cmd=['work/magma_wide_sample_sum','--gene-results',str(root/(trait+'.genes.raw')),'--gene-covar',str(matrix),'--model','condition-hide=Average','direction=greater','--settings','abbreviate=0','--out',str(prefix)]
 if not o.exists():
  p=subprocess.run(cmd,capture_output=True,text=True);pathlib.Path(str(prefix)+'.stdout').write_text(p.stdout+p.stderr);assert p.returncode==0 and o.exists(),(trait,p.returncode,p.stderr)
 with o.open() as f:lines=[x.split() for x in f if x.strip() and not x.startswith('#')]
 hdr=lines[0]
 for line in lines[1:]:
  x=dict(zip(hdr,line));p=None if x['P']=='NA' else float(x['P']);assert p is None or math.isfinite(p) and 0<=p<=1
  out.append({'trait':trait,'dataset':'HPA_human_skeletal_muscle_GSE233882','species':'human','cell_type_source_label':x['VARIABLE'],'n_genes':x['NGENES'],'beta':x['BETA'],'se':x['SE'],'p':p,'bh_q_100_trait_cell_tests':''})
 receipt.append({'trait':trait,'command':cmd,'result_sha256':hashlib.sha256(o.read_bytes()).hexdigest()})
assert len(out)==100
valid=sorted((x for x in out if x['p'] is not None),key=lambda x:x['p']);q=1
for i in range(len(valid),0,-1):q=min(q,valid[i-1]['p']*len(out)/i);valid[i-1]['bh_q_100_trait_cell_tests']=q
with open('outputs/stage4_hpa_muscle_four_trait_cell_properties.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
json.dump({'source':'https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data','selection':'post hoc human skeletal muscle sensitivity after first six matrices; 25 cluster tests across four trait models; BH family 100 planned tests; missing are retained in denominator','runs':receipt},open('outputs/stage4_hpa_muscle_gene_property_receipt.json','w'),indent=2)
print('rows',len(out),'estimable',len(valid),'BH<.05',sum(float(x['bh_q_100_trait_cell_tests'])<.05 for x in valid))
for x in sorted(valid,key=lambda x:x['p'])[:10]:print(x['trait'],x['cell_type_source_label'],x['p'],x['bh_q_100_trait_cell_tests'])
