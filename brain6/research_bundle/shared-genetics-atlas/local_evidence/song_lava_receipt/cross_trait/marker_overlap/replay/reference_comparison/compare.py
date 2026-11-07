import csv, hashlib, json, math
from pathlib import Path
import numpy as np
import pysam

p=Path(__file__).resolve().parent
panel=Path('/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel')
eur={r['sample'] for r in csv.DictReader(panel.open(),delimiter='\t') if r['super_pop']=='EUR'}
assert len(eur)==503
markers=list(csv.DictReader((p/'markers.tsv').open(),delimiter='\t')); assert len(markers)==60
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);samples=sorted(eur & set(vcf.header.samples));assert len(samples)==503
calls={};status=[]
for m in markers:
    records=[x for x in vcf.fetch('3',int(m['POS'])-1,int(m['POS'])) if x.pos==int(m['POS']) and m['SNP'] in (x.id or '').split(';') and m['A1'] in x.alleles and m['A2'] in x.alleles]
    if len(records)!=1:
        status.append({'SNP':m['SNP'],'status':'missing_or_ambiguous','matches':len(records)});continue
    x=records[0];a=x.alleles.index(m['A1']);b=x.alleles.index(m['A2']);arr=[]
    for sample in samples:
        gt=x.samples[sample].get('GT')
        arr.append(sum(z==a for z in gt) if gt and len(gt)==2 and all(z in (a,b) for z in gt) else np.nan)
    arr=np.asarray(arr,dtype=float);n=int(np.isfinite(arr).sum())
    if n<450 or np.nanstd(arr)==0:
        status.append({'SNP':m['SNP'],'status':'insufficient_or_invariant','matches':1,'n':n});continue
    calls[m['SNP']]=arr;freq=float(np.nanmean(arr)/2);ukb=float(m['FREQ'])
    status.append({'SNP':m['SNP'],'status':'matched','matches':1,'n':n,'a1_frequency':freq,'ukb_freq_field':ukb,
                   'frequency_difference_orientation_free':min(abs(freq-ukb),abs(1-freq-ukb))})
    print(m['SNP'],n,flush=True)
with (p/'ukb_ld.tsv').open() as f:
    ukb={row['SNP']:row for row in csv.DictReader(f,delimiter='\t')}
out=[]
for i,a in enumerate(markers):
    if a['SNP'] not in calls:continue
    for b in markers[i+1:]:
        if b['SNP'] not in calls:continue
        x,y=calls[a['SNP']],calls[b['SNP']];good=np.isfinite(x)&np.isfinite(y)
        if good.sum()<450 or np.std(x[good])==0 or np.std(y[good])==0:continue
        r=float(np.corrcoef(x[good],y[good])[0,1]);u=float(ukb[a['SNP']][b['SNP']]);dist=abs(int(a['POS'])-int(b['POS']))
        out.append({'snp1':a['SNP'],'snp2':b['SNP'],'distance_bp':dist,'n_eur':int(good.sum()),'ukb_r':u,'eur_r':r,
                    'ukb_r2':u*u,'eur_r2':r*r,'abs_r2_diff':abs(u*u-r*r)})
assert out
with (p/'pairwise.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
with (p/'marker_status.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['SNP','status','matches','n','a1_frequency','ukb_freq_field','frequency_difference_orientation_free'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(status)
summaries={}
for name,subset in [('all',out),('within_100kb',[x for x in out if x['distance_bp']<=100000]),('ukb_r2_ge_0.2',[x for x in out if x['ukb_r2']>=.2])]:
    if subset:
        d=[x['abs_r2_diff'] for x in subset];summaries[name]={'pairs':len(d),'median_abs_r2_diff':float(np.median(d)),'max_abs_r2_diff':max(d),'over_0.1':sum(x>.1 for x in d)}
receipt={'selected_markers':60,'matched_markers':len(calls),'individuals_in_EUR_panel':len(samples),'summaries':summaries,
 'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'url':url,'full_vcf_hashed':False,
 'note':'LD based on pairwise complete diploid dosages, r² avoids sign-orientation assumption; selected post-hoc, no genotypes saved.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
