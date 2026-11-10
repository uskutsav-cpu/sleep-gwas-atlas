import pandas as pd, json
from pathlib import Path
p=Path('work/meis1_gtex');o=Path('outputs')
d=pd.read_csv(p/'meis1_hfrspal_input.csv');counts={'input':len(d)}
for pop in ['nfe','fin']:
 x=pd.read_csv(o/f'stage6_meis1_hfrspal_{pop}_variant_index.csv')[['rs','gnomad_index','gnomad_af_alt','sign_to_anchor']].rename(columns={'rs':'id','gnomad_index':pop+'_idx','gnomad_af_alt':pop+'_af','sign_to_anchor':pop+'_sign'})
 d=d.merge(x,on='id',validate='one_to_one');counts[pop+'_merged']=len(d)
 d[pop+'_aligned_af']=d[pop+'_af'].where(d[pop+'_sign']==1,1-d[pop+'_af'])
d=d.merge(pd.read_csv(p/'meis1_hfrspal_frequencies.csv'),on='id',validate='one_to_one')
d['qtl_diff']=(d.qtl_alt_af-d.nfe_aligned_af).abs();d['hfrs_diff']=(d.hfrs_alt_af-d.fin_aligned_af).abs();counts['before_qc']=len(d)
d=d[(d.qtl_diff<=.2)&(d.hfrs_diff<=.2)].sort_values('pos38').reset_index(drop=True);counts['retained']=len(d);counts['added_palindromes']=d[d.rsid.isin(['rs12621789','rs6546230'])][['rsid','qtl_diff','hfrs_diff','nfe_sign','fin_sign']].to_dict('records');counts['max_qtl_diff']=d.qtl_diff.max();counts['max_hfrs_diff']=d.hfrs_diff.max()
d.to_csv(p/'meis1_hfrspal_aligned.csv',index=False);(o/'stage6_meis1_hfrspal_align_receipt.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts))
