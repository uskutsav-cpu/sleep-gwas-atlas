import csv,hashlib,json,collections,itertools,re,openpyxl,numpy as np
from pathlib import Path
book=Path('work/retallick_neural/published_supp2.xlsx'); out=Path('outputs'); w=openpyxl.load_workbook(book,read_only=True,data_only=True)
seq={};conflicting=set()
for name,sequence,length in w['SI Data 2.9'].iter_rows(min_row=2,values_only=True):
 if isinstance(name,str) and isinstance(sequence,str):
  if name in seq and seq[name]!=sequence.upper():conflicting.add(name)
  else:seq[name]=sequence.upper()
for name in conflicting:seq.pop(name,None)
s=w['SI_Data 2.12']; it=s.iter_rows(values_only=True);hdr=next(it);groups=collections.defaultdict(list)
for row in it:
 r=dict(zip(hdr,row)); ident=r['ID'];typ=r['type'];
 if isinstance(ident,str) and typ in ('ref','alt'):groups[ident].append(r)
kmers=[''.join(x) for n in (1,2,3) for x in itertools.product('ACGT',repeat=n)]
def kcounts(s):
 c=collections.Counter();
 for n in (1,2,3):c.update(s[i:i+n] for i in range(len(s)-n+1))
 return np.array([c[k] for k in kmers],float)
xs=[];ys=[];chrs=[];excluded=collections.Counter()
for ident,rs in groups.items():
 if len(rs)!=2 or {r['type'] for r in rs}!={'ref','alt'}:excluded['not_unique_pair']+=1;continue
 labels={r['emVar.group'] for r in rs}
 if len(labels)!=1:excluded['label_disagreement']+=1;continue
 label=next(iter(labels));
 if label not in ('fdr <= 0.05','n.s'):excluded['excluded_label']+=1;continue
 a,b=(next(r for r in rs if r['type']==t) for t in ('ref','alt'));chrmatch=re.fullmatch(r'chr(\d+):rs\d+:\d+',ident)
 if not chrmatch or not 1<=int(chrmatch[1])<=22:excluded['non_autosomal_id']+=1;continue
 sa=seq.get(a['CRS']);sb=seq.get(b['CRS']);
 if not sa or not sb or len(sa)!=199 or len(sb)!=199 or set(sa+sb)-set('ACGT'):excluded['sequence_unusable']+=1;continue
 if sum(x!=y for x,y in zip(sa,sb))!=1:excluded['not_single_base_pair']+=1;continue
 xs.append(np.abs(kcounts(sa)-kcounts(sb)));ys.append(int(label=='fdr <= 0.05'));chrs.append(int(chrmatch[1]))
X=np.stack(xs);y=np.array(ys);chrs=np.array(chrs);test=np.isin(chrs,[1,2]);tr=~test
assert test.any() and tr.any() and len(np.unique(y[test]))==2 and len(np.unique(y[tr]))==2
mean=X[tr].mean(0);sd=X[tr].std(0);sd[sd==0]=1;X=(X-mean)/sd;X=np.clip(X,-20,20)
# Features are fixed; only training chromosomes contribute to fitting.
wgt=np.zeros(X.shape[1]);intercept=0.;lam=.1
for _ in range(4000):
 z=np.clip(X[tr]@wgt+intercept,-30,30);err=1/(1+np.exp(-z))-y[tr];wgt-=.05*(X[tr].T@err/tr.sum()+lam*wgt);intercept-=.05*err.mean()
pred=1/(1+np.exp(-np.clip(X[test]@wgt+intercept,-30,30)));truth=y[test]; order=np.argsort(-pred,kind='stable');labs=truth[order];scores=pred[order];n_pos=int(sum(labs));n_neg=len(labs)-n_pos
# Group ties before computing AP and AUC; predictions are constant within each tie group.
ends=np.r_[np.flatnonzero(np.diff(scores)!=0),len(scores)-1]; tp=np.cumsum(labs)[ends];fp=(ends+1)-tp;rec=tp/n_pos;prec=tp/(tp+fp);ap=float(np.sum(np.diff(np.r_[0,rec])*prec));auc=float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/n_neg]))
res={'dataset':'Retallick-Townsley 2026 Supplementary Data 2 baseline and library','source_url':'https://doi.org/10.1038/s41467-026-77114-x','source_workbook_sha256':hashlib.sha256(book.read_bytes()).hexdigest(),'pair_group_count':len(groups),'excluded_reasons':dict(excluded),'train_n':int(tr.sum()),'train_detected':int(y[tr].sum()),'holdout_chromosomes':[1,2],'test_n':int(test.sum()),'test_detected':n_pos,'test_prevalence':float(n_pos/len(labs)),'test_auroc':auc,'test_average_precision':ap,'features':'absolute allele-pair differences in ACGT 1-3-mer counts; train-only standardization','model':'L2 logistic regression, lambda 0.1, 4000 fixed iterations, learning rate 0.05','target':'author FDR <=0.05 versus nonsignificant assay call; not a functional vs nonfunctional truth label','caveats':['Exploratory same-study holdout, no external donor or study validation.','Selected variants and binary statistical detection depend on assay power and selection; nonsignificant is not a biological negative.','This task is not disease prediction and cannot validate insomnia-frailty scores.','No AlphaGenome outputs or derivatives used.','Article CC BY-NC-ND: review reuse and redistribution terms.']}
(out/'stage6_neuronal_mpra_kmer_baseline.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
