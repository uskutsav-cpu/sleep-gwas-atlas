import csv,hashlib,json,math,collections
from pathlib import Path
p=Path(__file__).resolve().parent;source=Path('/Volumes/Extreme SSD/brain6-hfrs-candidate-results-20261028.tsv');source_receipt=Path('/Volumes/Extreme SSD/brain6-hfrs-candidate-receipt-20261028.json')
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(source)==json.loads(source_receipt.read_text())['output_sha256']
fi=list(csv.DictReader((p/'results.tsv').open(),delimiter='\t'));h=[r for r in csv.DictReader(source.open(),delimiter='\t') if r['score']=='original'];assert len(fi)==850 and len(h)==425 and len({r['rsid'] for r in h})==425
by={r['rsid']:r for r in h};candidate={}
for r in fi:candidate.setdefault(r['variant'],{})[r['trait']]=r
rows=[]
for snp,t in candidate.items():
 assert set(t)=={'FI','FFS'};a=t['FI'];b=t['FFS'];r=by.get(snp);rec={'variant':snp,'region':a['region'],'p_FI':a['p'],'p_FFS':b['p'],'status':'HFRS_MISSING'}
 if r:
  status=r['status'];coord='chr'+a['chr37']+':'+a['pos37']
  if status!='matched':rec['status']='HFRS_'+status
  elif r['grch37']!=coord or {r['a1'],r['a2']}!={a['effect_allele'],a['other_allele']}:rec['status']='HFRS_KEY_MISMATCH'
  else:
   try:
    vals=[float(a['p']),float(b['p']),float(r['p'])];assert all(math.isfinite(x) and 0<x<=1 for x in vals);mx=max(vals);rec.update(status='MATCHED',p_HFRS=vals[2],max_p=mx,max_p_bonferroni_425=min(1,425*mx),screen_hit=mx<.05/425)
   except (ValueError,AssertionError):rec['status']='INVALID_P'
 rows.append(rec)
keys=['variant','region','status','p_FI','p_FFS','p_HFRS','max_p','max_p_bonferroni_425','screen_hit']
with (p/'three_score.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
hits=[r for r in rows if r.get('screen_hit')];receipt={'hfrs_local_table_sha256':sha(source),'hfrs_original_receipt_sha256':sha(source_receipt),'fi_ffs_table_sha256':sha(p/'results.tsv'),'candidates':len(rows),'status_counts':dict(collections.Counter(r['status'] for r in rows)),'three_trait_maxp_hits':len(hits),'hit_regions':dict(collections.Counter(r['region'] for r in hits)),'hfrs_table_matches_original_receipt':True,'note':'Receipt binds the existing table, not an independently reconstructed raw source extraction.'};(p/'three_score_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
