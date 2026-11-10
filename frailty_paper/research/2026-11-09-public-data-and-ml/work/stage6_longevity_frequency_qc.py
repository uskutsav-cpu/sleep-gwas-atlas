import csv,json,numpy as np
from pathlib import Path
hits=json.loads(Path('outputs/stage6_longevity_candidate_screen.json').read_text())['matched']
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')));x=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True);assert list(x['variants'])==[r['variant'] for r in rows];f15={r['variant']:x['genotypes'][:,i].mean()/2 for i,r in enumerate(rows)}
f22={r['variant']:r for r in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv'))};out=[]
for h in hits:
 if h['score']=='15':freq=float(f15[h['variant']])
 else:
  v=f22[h['variant']];p=float(v['reference_EUR_frequency']);freq=p if v['counted_allele']==h['counted_allele'] else 1-p
 diff=abs(freq-h['eaf_counted']);out.append({'score':h['score'],'variant':h['variant'],'absolute_difference':diff,'reference_counted_frequency':freq,'target_counted_frequency':h['eaf_counted']})
res={s:{'n':sum(x['score']==s for x in out),'max_difference':max(x['absolute_difference'] for x in out if x['score']==s),'over_0_1':[x['variant'] for x in out if x['score']==s and x['absolute_difference']>0.1]} for s in ['15','22']}
Path('outputs/stage6_longevity_frequency_qc.json').write_text(json.dumps({'summary':res,'rows':out,'caveat':'Comparison is allele-coding QC only, not proof of strand, allele identity, or ancestral matching; target is meta-analysis with potentially heterogeneous cohorts and rounded frequencies.'},indent=2)+'\n');print(json.dumps(res,indent=2))
