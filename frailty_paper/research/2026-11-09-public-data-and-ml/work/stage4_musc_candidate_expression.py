import anndata as ad,numpy as np,pandas as pd,csv,collections,pathlib
q=list(csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv'))); genes=collections.defaultdict(set)
for r in q:
 g=r['gene ID'].split('.')[0]
 if g.startswith('ENSG'):genes[g].add(r['Candidate gene'])
a=ad.read_h5ad('work/SKM_MuSC_human_2023-06-22.h5ad',backed='r'); v=a.var.copy(); ids=v['ENSEMBL'].astype(str);counts=ids.value_counts(); eligible=[g for g in genes if counts.get(g,0)==1];missing=[g for g in genes if counts.get(g,0)==0];ambiguous=[g for g in genes if counts.get(g,0)>1];cols=np.array([int(np.flatnonzero(ids.to_numpy()==g)[0]) for g in eligible]); print('candidate genes',len(genes),'uniquely measured',len(eligible),'absent',len(missing),'ambiguous',len(ambiguous),flush=True)
from scipy.sparse import vstack
# Read small row blocks from the backed matrix and retain only candidate columns.
# This avoids loading the complete atlas expression matrix into memory.
X=vstack([a.X[i:i+500,:][:,cols] for i in range(0,a.n_obs,500)],format='csr');obs=a.obs;detailed=[]
for group in sorted(obs['annotation_level2'].unique()):
 for donor in sorted(obs['DonorID'].unique()):
  idx=np.flatnonzero((obs['annotation_level2']==group)&(obs['DonorID']==donor));n=len(idx)
  if n<10:continue
  sub=X[idx];means=np.asarray(sub.mean(axis=0)).ravel();det=np.asarray((sub>0).sum(axis=0)).ravel()/n
  for j,g in enumerate(eligible):detailed.append({'gene_id':g,'gene_symbols':'|'.join(sorted(genes[g])),'cell_group':group,'donor':donor,'n_cells':n,'mean_normalized_expression':float(means[j]),'fraction_detected':float(det[j])})
d=pd.DataFrame(detailed);out=[]
for (g,group),r in d.groupby(['gene_id','cell_group'],sort=True):
 out.append({'gene_id':g,'gene_symbols':r.gene_symbols.iloc[0],'cell_group':group,'n_donors_at_least_10_cells':len(r),'n_cells_in_those_donors':int(r.n_cells.sum()),'equal_donor_mean_normalized_expression':float(r.mean_normalized_expression.mean()),'median_donor_fraction_detected':float(r.fraction_detected.median()),'donors_with_detected_expression':int((r.fraction_detected>0).sum()),'evidence':'measured expression only; not genetic association or causality'})
pd.DataFrame(out).to_csv('outputs/stage4_human_musc_candidate_gene_expression.csv',index=False);pd.DataFrame(detailed).to_csv('work/stage4_human_musc_candidate_expression_by_donor.csv',index=False)
pathlib.Path('outputs/stage4_human_musc_gene_expression_coverage.json').write_text(__import__('json').dumps({'candidate_gene_ids':len(genes),'unique_ids_in_matrix':len(eligible),'absent_gene_ids':missing,'duplicated_gene_ids':ambiguous,'cell_subtypes':sorted(obs['annotation_level2'].unique().tolist()),'rule':'at least 10 cells per donor per subtype, equal donor weight; signal is normalized X, not genotype'},indent=2)+'\n')
print('rows',len(out),flush=True)
