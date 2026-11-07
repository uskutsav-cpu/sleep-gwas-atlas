import csv,json,math,hashlib
from pathlib import Path
import numpy as np,pysam
p=Path(__file__).resolve().parent
panel=Path('/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel')
eur={r['sample'] for r in csv.DictReader(panel.open(),delimiter='\t') if r['super_pop']=='EUR'};assert len(eur)==503
markers=list(csv.DictReader((p/'aligned_fi.tsv').open(),delimiter='\t'));assert len(markers)==2741
posmap={}
for m in markers:posmap.setdefault(int(m['POS']),[]).append(m)
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);samples=sorted(eur & set(vcf.header.samples));assert len(samples)==503;vcf.subset_samples(samples)
matches={};record_count=0
for x in vcf.fetch('3',min(posmap)-1,max(posmap)):
    if x.pos not in posmap:continue
    for m in posmap[x.pos]:
        if m['SNP'] not in (x.id or '').split(';') or m['A1'] not in x.alleles or m['A2'] not in x.alleles:continue
        a=x.alleles.index(m['A1']);b=x.alleles.index(m['A2']);vals=[]
        for sample in samples:
            gt=x.samples[sample].get('GT');vals.append(sum(z==a for z in gt) if gt and len(gt)==2 and all(z in (a,b) for z in gt) else np.nan)
        matches.setdefault(m['SNP'],[]).append(np.asarray(vals,dtype=float));record_count+=1
cols=[];selected=[];missing=[]
for m in markers:
    a=matches.get(m['SNP'],[])
    if len(a)!=1:missing.append({'SNP':m['SNP'],'reason':'missing_or_ambiguous','matches':len(a)});continue
    x=a[0];n=int(np.isfinite(x).sum())
    if n<450 or np.nanstd(x)==0:missing.append({'SNP':m['SNP'],'reason':'missing_calls_or_invariant','matches':1});continue
    x=np.where(np.isfinite(x),x,np.nanmean(x));x=(x-x.mean())/x.std(ddof=1)
    cols.append(x);selected.append(m)
assert len(selected)>3
X=np.column_stack(cols);del cols,matches
print('matched',len(selected),'excluded',len(missing),'starting SVD',flush=True)
U,d,Vt=np.linalg.svd(X,full_matrices=False);lam=d*d/(len(samples)-1)
keep=np.flatnonzero(np.cumsum(lam)/np.max(np.cumsum(lam))*100>=99)[0]+1
assert keep<503 and lam[keep-1]>0
z=np.array([float(m['STAT']) for m in selected]);ns=np.array([float(m['N']) for m in selected]);N=float(ns.mean())
corr=z/np.sqrt(z*z+ns-2);delta=(Vt[:keep,:]@corr)/np.sqrt(lam[:keep]);dtd=float(delta@delta)
sigma=(1-dtd)/(N-keep-1);omega=dtd/keep-sigma;h2=dtd-sigma*keep
with (p/'support.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
with (p/'excluded.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['SNP','reason','matches'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(missing)
out={'reference':'1000 Genomes Phase 3 EUR dosage, per-marker mean imputation','url':url,'whole_vcf_hashed':False,
 'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'n_people':len(samples),'requested_markers':len(markers),
 'matched_markers':len(selected),'excluded_markers':len(missing),'n_records_matched':record_count,
 'K_99':int(keep),'N_fi':N,'delta_sumsq':dtd,'sigma':sigma,'fi_omega':omega,'fi_h2_uncapped':h2,
 'min_retained_eigenvalue':float(lam[keep-1]),'support_sha256':hashlib.sha256((p/'support.tsv').read_bytes()).hexdigest(),
 'limitations':'Single-locus univariate sensitivity; reference small and non-study-matched; no individual genotypes saved.'}
(p/'eur_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
