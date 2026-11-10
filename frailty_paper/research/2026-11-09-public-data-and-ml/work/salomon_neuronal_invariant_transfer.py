import csv,gzip,json,collections,hashlib,runpy
from pathlib import Path
import numpy as np
source=Path('work/salomon_neuronal_variant_readout.tsv.gz')
# Select study chromosomes held out by the existing model. One source row per
# genomic contrast; exclude all duplicates instead of choosing an orientation.
groups=collections.defaultdict(list)
with gzip.open(source,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['chr'] in ('chr1','chr2'):groups[r['chr_pos_ref_alt']].append(r)
reasons=collections.Counter(); pairs=[]; group_labels=[]; frequency_categories=[]
for key,records in groups.items():
 if len(records)!=1:reasons['duplicate_contrast']+=1;continue
 r=records[0]; a=r['ref_seq'].upper();c=r['alt_seq'].upper(); label=r['NGN2_variant_is_significant']
 if r['NGN2_variant_has_readout']!='True' or label not in ('True','False'):reasons['missing_readout']+=1;continue
 if len(a)!=270 or len(c)!=270 or set(a+c)-set('ACGT') or sum(x!=y for x,y in zip(a,c))!=1:reasons['invalid_sequence_pair']+=1;continue
 ch,pos,ref,alt=key.split('-'); pos=int(pos); start=int(r['start']); end=int(r['end']); i=int(r['variant_pos']); strand=r['strand']
 if strand not in ('+','-') or end-start!=270:raise ValueError('Unexpected sequence coordinates: '+key)
 genomic_index=pos-1-start if strand=='+' else end-pos
 observed_index=next(j for j,(x,z) in enumerate(zip(a,c)) if x!=z)
 expected=(ref+alt) if strand=='+' else (ref+alt).translate(str.maketrans('ACGT','TGCA'))
 if i!=genomic_index or i!=observed_index or a[i]+c[i]!=expected:raise ValueError('Variant/sequence disagreement: '+key)
 pairs.append((a,c,int(label=='True')))
 group_labels.append(r['associated_gene_name'])
 frequency_categories.append(r['af_category'])
assert len(pairs)>1000 and len({p[2] for p in pairs})==2
# Reconstruct the fixed earlier training pipeline, without using assay data.
b=runpy.run_path('work/neuronal_mpra_kmer_baseline.py'); kmers=b['kmers']; idx={k:i for i,k in enumerate(kmers)};comp=str.maketrans('ACGT','TGCA');perm=np.array([idx[k.translate(comp)[::-1]] for k in kmers]);assert np.array_equal(perm[perm],np.arange(len(kmers)))
def invariant(x):return (x+x[...,perm])/2
X=invariant(np.stack(b['xs'])); y=np.array(b['ys']); train=~np.isin(np.array(b['chrs']),[1,2]);mean=X[train].mean(0);sd=X[train].std(0);sd[sd==0]=1;X=np.clip((X-mean)/sd,-20,20)
w=np.zeros(X.shape[1]);intercept=0.
for _ in range(4000):
 err=1/(1+np.exp(-np.clip(X[train]@w+intercept,-30,30)))-y[train];w-=.05*(X[train].T@err/train.sum()+.1*w);intercept-=.05*err.mean()
F=invariant(np.stack([np.abs(b['kcounts'](a)-b['kcounts'](c)) for a,c,_ in pairs]));labels=np.array([v for _,_,v in pairs]);scores=1/(1+np.exp(-np.clip(np.clip((F-mean)/sd,-20,20)@w+intercept,-30,30)))
# Tie-aware ROC and AP by score group.
order=np.argsort(-scores,kind='stable');s=scores[order];lab=labels[order];ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1];tp=np.cumsum(lab)[ends];fp=(ends+1)-tp;rec=tp/lab.sum();prec=tp/(ends+1);auc=float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/(len(lab)-lab.sum())]));ap=float(np.sum(np.diff(np.r_[0,rec])*prec))
receipt={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'training_workbook_sha256':b['res']['source_workbook_sha256'],'dataset':'Salomon NGN2 readouts chromosome 1 and 2; singleton source records by genomic contrast, exact 270-base one-SNV pairs','excluded':dict(reasons),'variant_and_strand_consistent_pairs':len(pairs),'train_n':int(train.sum()),'test_n':len(pairs),'test_significant':int(labels.sum()),'prevalence':float(labels.mean()),'auroc':auc,'average_precision':ap,'warning':'Different neuronal MPRA designs, sequence lengths and significance thresholds; no functional ground truth, disease or score outcomes. Previously redesigned baseline and posthoc benchmark; sequence/locus overlap and reuse conditions not fully assessed. No AlphaGenome predictions used.'}
Path('outputs/stage6_salomon_invariant_transfer.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
