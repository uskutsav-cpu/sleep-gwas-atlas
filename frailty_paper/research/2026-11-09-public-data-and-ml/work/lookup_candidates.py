import pathlib,csv,gzip,hashlib,json,sys,time
base=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1')
paths={'insomnia':base/'analysis-workspace/frailty_v1/harmonized/insomnia.harmonized.tsv.gz','FI':base/'analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz','FFS':base/'frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'}
source=list(csv.DictReader(open('outputs/published_insomnia_fi_124_cross_trait_candidates.csv')))
ids={x['variant_id'] for x in source}; assert len(ids)==124
found={}; metadata={}
for label,path in paths.items():
 start=time.time(); acc={}; count=0; opener=gzip.open if path.suffix=='.gz' else open
 with opener(path,'rt',newline='') as f:
  reader=csv.DictReader(f,delimiter='\t')
  for row in reader:
   count+=1
   if row['SNP'] in ids:acc.setdefault(row['SNP'],[]).append(row)
 found[label]=acc
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 metadata[label]={'path':str(path),'size_bytes':path.stat().st_size,'sha256':h.hexdigest(),'rows_read':count,'candidate_ids_found':len(acc),'ids_with_multiple_rows':sum(len(v)>1 for v in acc.values()),'seconds':round(time.time()-start,1)}
 print(label,metadata[label],flush=True)
result=[]
for x in source:
 out={'candidate':x['variant_id'],'song_chr':x['chromosome_reported'],'song_pos_reported':x['position_reported'],'song_allele':x['effect_allele_reported'],'song_other':x['other_allele_reported'],'song_rank_cpassoc':x['rank_by_cpassoc_p'],'song_cpassoc_p':x['cpassoc_p']}
 selected={}
 for label in paths:
  mm=found[label].get(x['variant_id'],[]);out[label+'_rows']=len(mm)
  if len(mm)==1:
   r=mm[0];selected[label]=r
   for field in ('CHR' if label!='FFS' else 'CHROM','BP' if label!='FFS' else 'POS','A1','A2','BETA','SE','P'):out[label+'_'+field]=r[field]
  else:
   for field in ('CHR' if label!='FFS' else 'CHROM','BP' if label!='FFS' else 'POS','A1','A2','BETA','SE','P'):out[label+'_'+field]=''
 if all(k in selected for k in ('insomnia','FI')):
  a,b=selected['insomnia'],selected['FI']; coord=(a['CHR']==b['CHR'] and a['BP']==b['BP']); pair={a['A1'],a['A2']}; pal=pair in ({'A','T'},{'C','G'})
  mode=('same' if (a['A1'],a['A2'])==(b['A1'],b['A2']) else 'swapped' if (a['A1'],a['A2'])==(b['A2'],b['A1']) else 'incompatible')
  out['FI_insomnia_coord_match']=str(coord);out['FI_insomnia_allele_relation']=mode;out['palindromic']=str(pal)
  if coord and not pal and mode in ('same','swapped'):
   fi=float(b['BETA'])*(1 if mode=='same' else -1); ins=float(a['BETA']); out['FI_beta_aligned_to_insomnia_A1']=fi;out['FI_insomnia_sign_relation']='same' if ins*fi>0 else 'opposite' if ins*fi<0 else 'zero'
  else:out['FI_beta_aligned_to_insomnia_A1']='';out['FI_insomnia_sign_relation']='not_interpreted'
 else:
  for k in ('FI_insomnia_coord_match','FI_insomnia_allele_relation','palindromic','FI_beta_aligned_to_insomnia_A1','FI_insomnia_sign_relation'):out[k]=''
 out['caveat']='Descriptive post-publication lookup, not independent replication or proof of shared causal variant; FFS A1 convention inferred; no risk score'
 result.append(out)
with open('outputs/published_candidate_marginal_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0].keys());w.writeheader();w.writerows(result)
pathlib.Path('work/candidate_lookup_receipt.json').write_text(json.dumps(metadata,indent=2))
print('done',len(result),flush=True)
