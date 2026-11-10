import runpy,numpy as np,json
from pathlib import Path
m=runpy.run_path('work/salomon_neuronal_invariant_transfer.py'); cats=np.asarray(m['frequency_categories']); y=m['labels'];p=m['scores'];results=[]
for cat in sorted(set(cats)):
 idx=cats==cat;labels=y[idx];scores=p[idx];pos=int(labels.sum());neg=len(labels)-pos
 row={'source_frequency_category':cat,'n':len(labels),'significant':pos,'prevalence':float(labels.mean())}
 if pos and neg:
  order=np.argsort(-scores,kind='stable');s=scores[order];v=labels[order];ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1];tp=np.cumsum(v)[ends];fp=ends+1-tp;rec=tp/pos;prec=tp/(ends+1)
  row['auroc']=float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/neg]));row['average_precision']=float(np.sum(np.diff(np.r_[0,rec])*prec))
 results.append(row)
receipt={'source_sha256':m['receipt']['source_sha256'],'results':results,'limitations':'Exploratory post-hoc subgroup analyses, no multiplicity correction or intervals; source frequency categories and labels not independently verified; selected neuronal assays, not outcomes or clinical validation.'}
Path('outputs/stage6_salomon_frequency_strata.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(results,indent=2))
