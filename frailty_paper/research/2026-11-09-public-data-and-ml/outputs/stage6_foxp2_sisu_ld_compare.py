#!/usr/bin/env python3
from pathlib import Path
import csv,gzip,json,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[1];w=root/'work';o=root/'outputs';data=pd.read_csv(w/'stage6_foxp2_susie_aligned.csv');idx=np.loadtxt(w/'stage6_foxp2_fin_reference_indices.txt',dtype=int)-1;g=np.fromfile(w/'stage6_foxp2_reference_genotypes_float64.bin',dtype='<f8').reshape(len(data),503);fin=g[:,idx];v=fin.var(axis=1);freq=fin.mean(axis=1)/2;mask=(v>0)&(freq>=.01)&(freq<=.99);data=data.loc[mask].reset_index(drop=True);g=g[mask];fin=fin[mask];bypos={int(r.pos):r for r in data.itertuples()};ref={};duplicates=set()
with gzip.open(w/'stage3_region_cache/FinnGen_insomnia.regional.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['#chrom']!='7':continue
  pos=int(r['pos']);s=bypos.get(pos)
  if s is None or set([r['ref'],r['alt']])!=set([s.anchor_effect_allele,s.anchor_other_allele]):continue
  key=f"7:{pos}:{r['ref']}:{r['alt']}"
  if s.rs in ref:duplicates.add(s.rs)
  else:ref[s.rs]=key
for rs in duplicates:ref.pop(rs,None)
rows=[];summary=[]
for anchor in ['rs35624276','rs12536335','rs62477310','rs2396682']:
 ai=np.where(data.rs==anchor)[0]
 if len(ai)!=1:raise ValueError(f'anchor absent {anchor}')
 ai=ai[0]
 def cor2(mat):
  a=mat[ai]-mat[ai].mean();b=mat-mat.mean(axis=1,keepdims=True);den=np.linalg.norm(b,axis=1)*np.linalg.norm(a);return ((b@a)/den)**2
 re=cor2(g);rf=cor2(fin);j=json.loads((o/f'stage6_foxp2_sisu_ld_{anchor}.json').read_text());sis={r['variation2']:r['r2'] for r in j['ld']};r=[]
 for k,s in enumerate(data.itertuples()):
  key=ref.get(s.rs)
  if key not in sis:continue
  z={'anchor':anchor,'variant':s.rs,'sisu_key':key,'sisu_r2':float(sis[key]),'eur503_r2':float(re[k]),'fin99_r2':float(rf[k])};r.append(z);rows.append(z)
 x=pd.DataFrame(r)
 summary.append({'anchor':anchor,'n_comparable':len(x),'eur_median_abs_r2_error':float((x.eur503_r2-x.sisu_r2).abs().median()),'fin99_median_abs_r2_error':float((x.fin99_r2-x.sisu_r2).abs().median()),'eur_mean_abs_r2_error':float((x.eur503_r2-x.sisu_r2).abs().mean()),'fin99_mean_abs_r2_error':float((x.fin99_r2-x.sisu_r2).abs().mean()),'eur_r2_correlation':float(x.eur503_r2.corr(x.sisu_r2)),'fin99_r2_correlation':float(x.fin99_r2.corr(x.sisu_r2))})
pd.DataFrame(rows).to_csv(o/'stage6_foxp2_sisu_anchor_r2_comparison.csv.gz',index=False,compression='gzip');pd.DataFrame(summary).to_csv(o/'stage6_foxp2_sisu_anchor_r2_summary.csv',index=False);print(pd.DataFrame(summary).to_string(index=False))
