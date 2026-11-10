import csv,hashlib,pathlib
root=pathlib.Path('outputs')
def read(s): return list(csv.DictReader((root/s).open()))
top=read('published_insomnia_fi_top30_candidates.csv'); lookup={r['candidate']:r for r in read('published_candidate_marginal_lookup.csv')}; proto={r['variant']:r for r in read('unvalidated_15_variant_score_prototype.csv')}; joint={r['variant']:r for r in read('mvp_insomnia_hfrs_joint_candidate_screen.csv')}
assert len(top)==30 and len(proto)==15 and set(proto)<=set(x['variant_id'] for x in top)
rows=[]
for x in top:
 v=x['variant_id']; y=lookup[v]; retained=v in proto
 rows.append({'song_rank':x['rank_by_cpassoc_p'],'variant':v,'song_chr':x['chromosome_reported'],'song_position_reported':x['position_reported'],'song_cpassoc_p':x['cpassoc_p'],'mhc_interval_flag':x['approx_chr6_MHC_flag'],'prototype_retained':retained,'insomnia_source_rows':y['insomnia_rows'],'fi_source_rows':y['FI_rows'],'ffs_source_rows':y['FFS_rows'],'FI_FFS_insomnia_allele_match_relation':y['FI_insomnia_allele_relation']+'/'+y['FFS_insomnia_allele_relation'] if retained else '', 'MVP_insomnia_HFRS_joint_p_exploratory':joint[v]['joint_p_max'] if retained else ''})
assert sum(x['prototype_retained'] for x in rows)==15
assert all((x['mhc_interval_flag']=='no') == x['prototype_retained'] for x in rows)
assert all(x['insomnia_source_rows']=='0' and x['fi_source_rows']=='0' for x in rows if not x['prototype_retained'])
p=root/'candidate_top30_selection_audit.csv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for x in ['published_insomnia_fi_top30_candidates.csv','published_candidate_marginal_lookup.csv','unvalidated_15_variant_score_prototype.csv','mvp_insomnia_hfrs_joint_candidate_screen.csv','candidate_top30_selection_audit.csv']: print(x,hashlib.sha256((root/x).read_bytes()).hexdigest())
