import csv, pathlib
p=pathlib.Path('outputs/published_candidate_marginal_lookup.csv')
with p.open(newline='') as f:rows=list(csv.DictReader(f))
for x in rows:
 for t,chrname,posname in [('insomnia','CHR','BP'),('FI','CHR','BP'),('FFS','CHROM','POS')]:
  x[t+'_same_reported_position_as_song']=str(x[t+'_rows']=='1' and x[t+'_'+chrname]==x['song_chr'] and x[t+'_'+posname]==x['song_pos_reported']) if x[t+'_rows']=='1' else ''
 x['FFS_insomnia_allele_relation']='';x['FFS_beta_aligned_to_insomnia_A1']='';x['FFS_insomnia_sign_relation']=''
 if x['insomnia_rows']=='1' and x['FFS_rows']=='1':
  a=(x['insomnia_A1'],x['insomnia_A2']); b=(x['FFS_A1'],x['FFS_A2']); relation='same' if a==b else 'swapped' if a==b[::-1] else 'incompatible';x['FFS_insomnia_allele_relation']=relation
  coord=x['insomnia_CHR']==x['FFS_CHROM'] and x['insomnia_BP']==x['FFS_POS'];pal=set(a) in ({'A','T'},{'C','G'})
  if coord and not pal and relation in ('same','swapped'):
   v=float(x['FFS_BETA'])*(1 if relation=='same' else -1); x['FFS_beta_aligned_to_insomnia_A1']=v
   product=v*float(x['insomnia_BETA']);x['FFS_insomnia_sign_relation']='same' if product>0 else 'opposite' if product<0 else 'zero'
  else:x['FFS_insomnia_sign_relation']='not_interpreted'
with p.open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
