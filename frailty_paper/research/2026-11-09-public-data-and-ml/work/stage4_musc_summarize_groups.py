import csv,pathlib
p=pathlib.Path('work/stage4_musc_scdrs_downstream/PASS_Insomnia_Jansen2019.scdrs_group.annotation_level2')
a=list(csv.DictReader(p.open(),delimiter='\t'));assert len(a)==4
for key,out in [('assoc_mcp','assoc_bh_fdr_4_subtypes'),('hetero_mcp','hetero_bh_fdr_4_subtypes')]:
 order=sorted(range(len(a)),key=lambda i:float(a[i][key]));q=1
 for rank in range(len(order),0,-1):
  i=order[rank-1];q=min(q,float(a[i][key])*len(order)/rank);a[i][out]=q
with open('outputs/stage4_human_musc_insomnia_scdrs_groups.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(a[0]));w.writeheader();w.writerows(a)
for x in sorted(a,key=lambda r:float(r['assoc_mcp'])):print(x['group'],x['n_cell'],x['assoc_mcp'],x['assoc_bh_fdr_4_subtypes'],x['hetero_bh_fdr_4_subtypes'])
