import csv,hashlib,json,pathlib,subprocess,math,re
root=pathlib.Path('work/stage4_magma');binary='work/magma_wide_sample_sum'
traits=['mak_fi','ye_ffs','mak_fi_brown','ye_ffs_brown','jansen_insomnia_brown','atkins_fi_brown']; datasets={'GSE67835_Human_Cortex_woFetal':('brain','human'),'GSE87544_Mouse_Hypothalamus':('brain','mouse mapped to human Ensembl IDs'),'PBMC_10x_68k':('immune','human'),'GSE98816_Mouse_Brain_Vascular':('vascular','mouse mapped to human Ensembl IDs'),'MouseCellAtlas_Muscle':('muscle','mouse mapped to human Ensembl IDs'),'TabulaMuris_FACS_Limb_Muscle':('muscle','mouse mapped to human Ensembl IDs')}
allrows=[];runs=[]
for trait in traits:
 raw=root/(trait+'.genes.raw')
 log=root/(trait+'.log'); geneout=root/(trait+'.genes.out')
 if not raw.is_file() or not log.is_file() or not geneout.is_file() or 'End time is' not in log.read_text():continue
 expected=re.findall(r'found (\d+) genes containing valid SNPs',log.read_text());assert len(expected)==1 and sum(1 for _ in geneout.open())-1==int(expected[0])
 for dataset,(domain,species) in datasets.items():
  matrix=root/'matrices'/(dataset+'.txt');prefix=root/(trait+'_'+dataset);o=pathlib.Path(str(prefix)+'.gsa.out')
  cmd=[binary,'--gene-results',str(raw),'--gene-covar',str(matrix),'--model','condition-hide=Average','direction=greater','--settings','abbreviate=0','--out',str(prefix)]
  if not o.exists():
   r=subprocess.run(cmd,text=True,capture_output=True);(root/(trait+'_'+dataset+'.stdout')).write_text(r.stdout+r.stderr)
   if r.returncode!=0 or not o.exists():raise RuntimeError((dataset,r.returncode,r.stdout[-1000:],r.stderr[-1000:]))
  with o.open() as f:
   lines=[x.split() for x in f if x.strip() and not x.startswith('#')]
  header=lines[0]; assert {'VARIABLE','NGENES','BETA','SE','P'}.issubset(header)
  for l in lines[1:]:
   x=dict(zip(header,l));p=None if x['P']=='NA' else float(x['P']);assert p is None or (math.isfinite(p) and 0<=p<=1)
   allrows.append({'trait':trait,'dataset':dataset,'domain':domain,'species':species,'cell_type':x['VARIABLE'],'n_genes':x['NGENES'],'beta':x['BETA'],'se':x['SE'],'p':p})
  runs.append({'trait':trait,'dataset':dataset,'command':cmd,'gsa_sha256':hashlib.sha256(o.read_bytes()).hexdigest()})
for trait in traits:
 rows=[r for r in allrows if r['trait']==trait];o=sorted((r for r in rows if r['p'] is not None),key=lambda x:x['p']);q=1
 # Preserve all planned hypotheses in the BH denominator; failed regressions have no p/q.
 for j in range(len(o),0,-1):q=min(q,o[j-1]['p']*len(rows)/j);o[j-1]['bh_q_all_six_matrices']=q
 for r in rows:
  if r['p'] is None:r['bh_q_all_six_matrices']=''
if allrows:
 with open('outputs/stage4_magma_cell_type_exploratory.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=allrows[0].keys());w.writeheader();w.writerows(allrows)
 pathlib.Path('outputs/stage4_magma_cell_type_receipt.json').write_text(json.dumps({'selection':'six exploratory matrices; all cell-type variables per trait and gene model in one BH family; model sensitivity and matrix choices not prespecified before analysis','runs':runs},indent=2)+'\n')
print('rows',len(allrows),'runs',len(runs))
