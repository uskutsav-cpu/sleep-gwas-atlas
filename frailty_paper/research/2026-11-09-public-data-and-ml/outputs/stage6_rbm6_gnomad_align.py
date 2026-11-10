from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('outputs');W=Path('work');d=pd.read_csv(W/'stage6_rbm6_gnomad_pair_input.csv');n=pd.read_csv(O/'stage6_rbm6_gnomad_nfe_variant_index.csv');f=pd.read_csv(O/'stage6_rbm6_gnomad_fin_variant_index.csv');d=d.merge(n[['rs','gnomad_index','gnomad_af_alt','sign_to_anchor']].rename(columns={'gnomad_index':'nfe_idx','gnomad_af_alt':'nfe_ref_af','sign_to_anchor':'nfe_sign'}),on='rs',validate='one_to_one').merge(f[['rs','gnomad_index','gnomad_af_alt','sign_to_anchor']].rename(columns={'gnomad_index':'fin_idx','gnomad_af_alt':'fin_ref_af','sign_to_anchor':'fin_sign'}),on='rs',validate='one_to_one');d['nfe_anchor_af']=np.where(d.nfe_sign==1,d.nfe_ref_af,1-d.nfe_ref_af);d['fin_anchor_af']=np.where(d.fin_sign==1,d.fin_ref_af,1-d.fin_ref_af)
j=pd.read_csv(W/'stage3_region_cache/Jansen_insomnia.regional.tsv.gz',sep='\t');j=j[j.CHR.astype(str).eq('3')];j['jkey']=j.SNP.astype(str)+'|'+j.BP.astype(str)+'|'+j.A1.astype(str)+'|'+j.A2.astype(str);assert j.jkey.is_unique
h=pd.read_csv(W/'stage6_rbm6_hfrs_2mb.tsv.gz',sep='\t').dropna(subset=['rs_id']);h['hkey']=h.rs_id.astype(str)+'|'+h.base_pair_location.astype(str)+'|'+h.effect_allele.astype(str)+'|'+h.other_allele.astype(str);assert h.hkey.is_unique
jd={r.jkey:float(r.FRQ) for r in j.itertuples()};hd={r.hkey:float(r.effect_allele_frequency) for r in h.itertuples()};jf=[];hf=[]
for r in d.itertuples():
 def lookup(D,keyprefix,a,b):
  k=keyprefix+'|'+a+'|'+b
  if k in D:return D[k]
  k=keyprefix+'|'+b+'|'+a
  if k in D:return 1-D[k]
  return np.nan
 jf.append(lookup(jd,f'{r.rs}|{int(r.pos37)}',r.anchor_effect_allele,r.anchor_other_allele));hf.append(lookup(hd,f'{r.rs}|{r.pos}',r.anchor_effect_allele,r.anchor_other_allele))
d['jansen_anchor_af']=jf;d['hfrs_anchor_af']=hf;d['jansen_af_diff']=abs(d.jansen_anchor_af-d.nfe_anchor_af);d['hfrs_af_diff']=abs(d.hfrs_anchor_af-d.fin_anchor_af);pre=len(d);d=d[d.jansen_af_diff.le(.2)&d.hfrs_af_diff.le(.2)].copy();d.to_csv(O/'stage6_rbm6_gnomad_aligned.csv',index=False);rec={'input':pre,'retained':len(d),'excluded_af_or_missing':pre-len(d),'mean_jansen_nfe_diff':float(d.jansen_af_diff.mean()),'mean_hfrs_fin_diff':float(d.hfrs_af_diff.mean()),'max_jansen_nfe_diff':float(d.jansen_af_diff.max()),'max_hfrs_fin_diff':float(d.hfrs_af_diff.max()),'contains_hfrs_lead':bool(d.rs.eq('rs2282749').any()),'jansen_original_lead_excluded_by_window':True};(O/'stage6_rbm6_gnomad_align_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec)
