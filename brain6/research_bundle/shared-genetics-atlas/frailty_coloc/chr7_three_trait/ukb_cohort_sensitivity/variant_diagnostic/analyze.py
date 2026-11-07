import csv, math, json, hashlib
from pathlib import Path
p=Path(__file__).resolve().parent; f=Path('/Volumes/Brain6Work/chr7-three-trait-ukb-sensitivity.tsv');r=list(csv.DictReader(f.open(),delimiter='\t'));assert len(r)==1055

def scores(pre,sd):
 x=[]
 for v in r:
  b=float(v[pre+'_beta']);s=float(v[pre+'_se']);w=sd**2/(sd**2+s*s);x.append(.5*(math.log1p(-w)+w*(b/s)**2))
 m=max(x);e=[math.exp(t-m) for t in x];ss=math.fsum(e);return [a/ss for a in e]
def credible(q):
 order=sorted(range(len(q)),key=lambda i:q[i],reverse=True);s=0;z=[]
 for i in order:
  z.append(i);s+=q[i]
  if s>=.95:break
 return set(z)
q={a:scores(a,sd) for a,sd in [('fi',.15),('ukb_fi',.15),('adhd',.2),('insomnia',.2)]}
cs={a:credible(v) for a,v in q.items()};result={'n':len(r),'input_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'trait_95_credible_set_sizes':{k:len(v) for k,v in cs.items()},'trait_top_variant':{k:r[max(range(len(v)),key=lambda i:v[i])]['rsid'] for k,v in q.items()},'credible_set_intersection_counts':{},'pair_comparison':{}}
for source in ['fi','ukb_fi']:
 for t in ['adhd','insomnia']:
  qa=q[source];qb=q[t];dot=math.fsum(a*b for a,b in zip(qa,qb));weights=[a*b/dot for a,b in zip(qa,qb)];best=max(range(len(r)),key=lambda i:weights[i]);c=credible(weights)
  result['credible_set_intersection_counts'][source+'_'+t]=len(cs[source]&cs[t]);result['pair_comparison'][source+'_'+t]={'normalized_bayes_factor_overlap':dot,'h4_h3_odds_at_p12_1e-5':1e3*dot/(1-dot),'h4_conditional_top_rsid':r[best]['rsid'],'h4_conditional_top_weight':weights[best],'h4_conditional_95_set_size':len(c)}

with (p.parent/'pairwise.tsv').open() as h:
 for row in csv.DictReader(h,delimiter='\t'):
  if float(row['p12'])!=1e-5: continue
  key=('fi' if row['source']=='meta' else 'ukb_fi')+'_'+row['pair'].lower()
  assert abs(result['pair_comparison'][key]['h4_h3_odds_at_p12_1e-5']-float(row['H4'])/float(row['H3']))<1e-10
assert all(abs(math.fsum(v)-1)<1e-12 for v in q.values())
result['h4_h3_odds_match_parent']=True
result['fi_source_95_set_intersection_count']=len(cs['fi']&cs['ukb_fi']);result['limitations']='All summaries are conditional on one causal signal per trait and selected support. No LD-based conditional test, replication, calibration or probability of gene causality.'
(p/'diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
