import collections,csv,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;r=list(csv.DictReader((p/'results.tsv').open(),delimiter='\t'));t=list(csv.DictReader((p/'three_score.tsv').open(),delimiter='\t'));assert len(r)==850 and len(t)==425
assert len({(x['variant'],x['trait']) for x in r})==850 and len({x['variant'] for x in t})==425
by={}
for x in r:
 assert x['status']=='MATCHED' and int(x['exact_rows'])==1 and int(x['same_id_rows'])>=1
 v=float(x['p']);assert 0<v<=1 and abs(float(x['p_bonferroni_850'])-min(1,850*v))<1e-12
 assert (x['screen_hit']=='True')==(v<.05/850);by.setdefault(x['variant'],{})[x['trait']]=x
A={s for s,x in by.items() if x['FI']['screen_hit']=='True'};B={s for s,x in by.items() if x['FFS']['screen_hit']=='True'};assert len(A)==129 and len(B)==71 and len(A&B)==71
for x in t:
 assert x['status']=='MATCHED';values=[float(x[k]) for k in ('p_FI','p_FFS','p_HFRS')];assert abs(float(x['max_p'])-max(values))<1e-14;assert abs(float(x['max_p_bonferroni_425'])-min(1,425*max(values)))<1e-12;assert (x['screen_hit']=='True')==(max(values)<.05/425)
assert not any(x['screen_hit']=='True' for x in t)
res={'n_rows':len(r),'n_variants':len(by),'FI_hits':len(A),'FFS_hits':len(B),'both_hits':len(A&B),'FI_only_threshold_pattern':len(A-B),'FFS_only_threshold_pattern':len(B-A),'three_trait_hits':0,'warning':'Arithmetic and keys verified; statistical calibration, selection and source origin not validated.'};(p/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
