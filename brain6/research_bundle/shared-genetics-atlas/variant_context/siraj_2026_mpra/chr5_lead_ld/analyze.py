import csv,collections,json,math,random,hashlib
from pathlib import Path
import pysam
out=Path(__file__).resolve().parent
panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
with panel.open() as f: eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
assert len(eur)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
variants=[('rs12658032',103904226,'A','G'),('rs2431108',103947968,'T','C'),('rs77960',103964585,'G','A')]
v=pysam.VariantFile(url);calls={};meta={}
for rs,pos,ref,alt in variants:
 hits=[r for r in v.fetch('5',pos-1,pos) if r.pos==pos and r.id and rs in r.id.split(';') and r.ref==ref and alt in r.alts]
 assert len(hits)==1,(rs,len(hits));r=hits[0];ai=r.alleles.index(alt);d={}
 for sample in sorted(eur & set(v.header.samples)):
  s=r.samples[sample];gt=s.get('GT')
  if gt and len(gt)==2 and s.phased and all(x in (0,ai) for x in gt):d[sample]=tuple(int(x==ai) for x in gt)
 calls[rs]=d;meta[rs]={'position37':pos,'ref':ref,'alt':alt,'n_people':len(d),'frequency_alt':sum(sum(x) for x in d.values())/(2*len(d))}
def correlation(pairs):
 n=2*len(pairs);a=sum(x[0] for x in pairs)/n;b=sum(x[1] for x in pairs)/n;ab=sum(x[2] for x in pairs)/n;return (ab-a*b)/math.sqrt(a*(1-a)*b*(1-b))
rows=[]
for i,(a,*_) in enumerate(variants):
 for b,*_ in variants[i+1:]:
  people=sorted(set(calls[a])&set(calls[b]));pairs=[(sum(calls[a][s]),sum(calls[b][s]),sum(x*y for x,y in zip(calls[a][s],calls[b][s]))) for s in people]
  r=correlation(pairs);rng=random.Random(20261027);bs=sorted(correlation([pairs[rng.randrange(len(pairs))] for _ in pairs])**2 for i in range(2000))
  rows.append({'variant_1':a,'alt_1':meta[a]['alt'],'variant_2':b,'alt_2':meta[b]['alt'],'people':len(people),'r':r,'r2':r*r,'r2_boot_025':bs[49],'r2_boot_975':bs[1949]})
with (out/'pairwise_ld.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'source_url':url,'whole_remote_vcf_sha256':'not_fetched','reference_panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'variants':meta,'bootstrap_replicates':2000,'seed':20261027,'output_sha256':hashlib.sha256((out/'pairwise_ld.tsv').read_bytes()).hexdigest()}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'variants':meta,'rows':rows},indent=2))
