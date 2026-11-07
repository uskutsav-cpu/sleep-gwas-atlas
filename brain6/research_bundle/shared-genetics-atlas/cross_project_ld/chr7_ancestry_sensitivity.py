"""Descriptive reference-population sensitivity, not a cohort-matched analysis."""
import csv, pathlib, json, hashlib, collections, math, pysam
p=pathlib.Path(__file__).resolve().parent
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
source=p.parent/'frailty_coloc/conditional_variant_posteriors.tsv'
with panel.open() as f: panel_rows=list(csv.DictReader(f,delimiter='\t'))
groups={k:{x['sample'] for x in panel_rows if (x['super_pop']=='EUR' if k=='EUR' else x['pop']==k)} for k in ('EUR','FIN','GBR','CEU','IBS','TSI')}
assert len(groups['EUR'])==503 and sum(len(groups[k]) for k in ('FIN','GBR','CEU','IBS','TSI'))==503
variants={'rs10953765':(114291435,'A','G'),'rs1456031':(114296102,'T','C'),'rs1476535':(114071035,'C','T'),'rs12536335':(114043159,'G','A')}
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);calls={}
for rs,(pos,counted,other) in variants.items():
    matches=[x for x in vcf.fetch('7',pos-1,pos) if x.pos==pos and x.id==rs and counted in x.alleles and other in x.alleles]
    assert len(matches)==1,(rs,len(matches))
    rec=matches[0]; a=rec.alleles.index(counted);b=rec.alleles.index(other);calls[rs]={}
    for sample in groups['EUR'] & set(vcf.header.samples):
        c=rec.samples[sample];gt=c.get('GT')
        if gt and len(gt)==2 and c.phased and all(z in (a,b) for z in gt): calls[rs][sample]=tuple(int(z==a) for z in gt)
freq=[];pairs=[]
for group,samples in groups.items():
    for rs,(pos,counted,other) in variants.items():
        available=samples & calls[rs].keys(); n=2*len(available);assert n>0
        freq.append({'group':group,'variant':rs,'pos_grch37':pos,'counted_allele':counted,'other_allele':other,'individuals':len(available),'allele_frequency':sum(sum(calls[rs][x]) for x in available)/n})
    for a,b in [('rs10953765','rs1456031'),('rs10953765','rs1476535'),('rs10953765','rs12536335'),('rs1456031','rs1476535'),('rs1456031','rs12536335')]:
        shared=sorted(samples & calls[a].keys() & calls[b].keys());n=2*len(shared)
        xs=[z for s in shared for z in calls[a][s]];ys=[z for s in shared for z in calls[b][s]]
        f1=sum(xs)/n;f2=sum(ys)/n;f11=sum(x*y for x,y in zip(xs,ys))/n
        den=math.sqrt(f1*(1-f1)*f2*(1-f2));assert den>0
        r=(f11-f1*f2)/den;pairs.append({'group':group,'variant_1':a,'variant_2':b,'individuals':len(shared),'r':r,'r2':r*r})
for name,data in [('chr7_ancestry_frequencies.tsv',freq),('chr7_ancestry_ld.tsv',pairs)]:
    with (p/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
(p/'chr7_ancestry_receipt.json').write_text(json.dumps({'reference_url':url,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'source_posterior_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'population_counts':{k:len(v) for k,v in groups.items()},'variants':variants,'method':'phased counted-allele frequency and haplotype correlation; 1000G subpopulations; no genotype data saved'},indent=2)+'\n')
print('Saved',len(freq),'frequency and',len(pairs),'pairwise group summaries')
