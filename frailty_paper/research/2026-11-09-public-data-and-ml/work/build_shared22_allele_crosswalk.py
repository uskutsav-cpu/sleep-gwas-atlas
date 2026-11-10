import csv,json,hashlib,pathlib
base=pathlib.Path('outputs'); p=base/'expanded_22_candidate_risk_gene_table.csv'; q=base/'shared22_eight_source_evidence_matrix.csv'
x=list(csv.DictReader(p.open())); y={r['variant']:r for r in csv.DictReader(q.open())}; assert len(x)==len(y)==22
out=[]
for r in x:
 a=y[r['variant']]; assert {r['counted_allele'],r['other_allele']}=={a['score_counted_allele'],a['other_allele']}; assert r['nearest_gene_body_positional']==a['nearest_gene_body_positional']
 flip=r['counted_allele']!=a['score_counted_allele']; sign=-1 if flip else 1
 for f,g in [('FI_beta_aligned','atkins_FI_beta_aligned'),('FFS_beta_aligned','ye_FFS_beta_aligned'),('HFRS_beta_aligned','hfrs_beta_aligned')]:
  assert abs(float(r[f])*sign-float(a[g]))<1e-12
 out.append({'variant':r['variant'],'chr_GRCh37_assumed':r['GRCh37_chr_assumed'],'position_GRCh37_assumed':r['GRCh37_pos_assumed'],'nearest_gene_body_positional':r['nearest_gene_body_positional'],'original_table_counted_allele':r['counted_allele'],'original_table_other_allele':r['other_allele'],'shared_count_counted_allele':a['score_counted_allele'],'shared_count_other_allele':a['other_allele'],'orientation_switched':str(flip),'FI_beta_shared_count_allele':a['atkins_FI_beta_aligned'],'FFS_beta_shared_count_allele':a['ye_FFS_beta_aligned'],'HFRS_beta_shared_count_allele':a['hfrs_beta_aligned'],'MVP_log_OR_shared_count_allele_approx':a['mvp_insomnia_beta_aligned']})
z=base/'shared22_allele_orientation_crosswalk.csv'
with z.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
print(json.dumps({'total':len(out),'switched':sum(r['orientation_switched']=='True' for r in out),'switched_ids':[r['variant'] for r in out if r['orientation_switched']=='True'],'input_hashes':[h(p),h(q)],'output_hash':h(z)},indent=2))
