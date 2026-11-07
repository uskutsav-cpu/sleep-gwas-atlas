import csv,json,math,hashlib,pathlib,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');chain=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
eur={r['sample'] for r in csv.DictReader(panel.open(),delimiter='\t') if r['super_pop']=='EUR'};assert len(eur)==503
rows=[r for r in csv.DictReader((p/'conditional_snps.tsv').open(),delimiter='\t') if r['credible95']=='TRUE'];vs={r['variant']:r['rsid'] for r in rows};lo=LiftOver(str(chain));url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);calls={};records={}
for var,rs in sorted(vs.items()):
 _,pos,ref,alt=var.split('_');converted=lo.convert_coordinate('chr5',int(pos)-1)
 if len(converted)!=1 or converted[0][0]!='chr5' or converted[0][2]!='+': records[rs]={'status':'lift_failed'};continue
 pos37=converted[0][1]+1;hits=[x for x in vcf.fetch('5',pos37-1,pos37) if x.pos==pos37 and x.id==rs and ref in x.alleles and alt in x.alleles]
 if len(hits)!=1:records[rs]={'status':'reference_match_not_unique','position37':pos37,'n':len(hits)};continue
 rec=hits[0];ia=rec.alleles.index(alt);ir=rec.alleles.index(ref);data={}
 for sample in sorted(eur & set(vcf.header.samples)):
  c=rec.samples[sample];g=c.get('GT')
  if g and len(g)==2 and all(t is not None and t in (ia,ir) for t in g) and c.phased:data[sample]=tuple(int(t==ia) for t in g)
 calls[rs]=data;records[rs]={'status':'matched','position37':pos37,'position38':int(pos),'ref':ref,'alt':alt,'n_individuals':len(data),'reference_alleles':list(rec.alleles)}
results=[];keys=sorted(calls)
for k,one in enumerate(keys):
 for two in keys[k+1:]:
  common=set(calls[one])&set(calls[two]);x=[];y=[]
  for sample in sorted(common):x.extend(calls[one][sample]);y.extend(calls[two][sample])
  n=len(x);a=sum(x)/n;b=sum(y)/n;ab=sum(c*d for c,d in zip(x,y))/n;den=math.sqrt(a*(1-a)*b*(1-b));r=(ab-a*b)/den if den else None
  results.append({'one':one,'two':two,'n_individuals':len(common),'r':r,'r2':r*r if r is not None else None})
with (p/'pairwise_ld.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['one','two','n_individuals','r','r2'],delimiter='\t');w.writeheader();w.writerows(results)
receipt={'reference':url,'selected_variant_records':records,'eur_n':len(eur),'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'chain_sha256':hashlib.sha256(chain.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((p/'conditional_snps.tsv').read_bytes()).hexdigest(),'note':'Remote indexed reference files not whole-file hashed; phased pairwise complete r²; not study-matched'}
(p/'reference_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print(json.dumps(results,indent=2))
