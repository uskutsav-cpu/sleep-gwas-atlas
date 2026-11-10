import csv,json,hashlib,math,pysam,numpy as np
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));v=pysam.VariantFile(url);vs={'rs7734045':(135538341,'C','A'),'rs7731706':(135570987,'G','T')};rec={}
for rs,(pos,a,b) in vs.items():
 x=[x for x in v.fetch('5',pos-1,pos) if x.id==rs and set(x.alleles)=={a,b} and len(x.alleles)==2];assert len(x)==1,(rs,len(x));rec[rs]=x[0]
o={}
for pop in ['AFR','EUR']:
 people=[r['sample'] for r in panel if r['super_pop']==pop and r['sample'] in v.header.samples];X=[];Y=[]
 for person in people:
  a,b=[rec[k].samples[person] for k in vs];g=a.get('GT');h=b.get('GT')
  if not (a.phased and b.phased and g and h and len(g)==2 and len(h)==2 and None not in g and None not in h):continue
  X.extend([int(z==rec['rs7734045'].alleles.index('C')) for z in g]);Y.extend([int(z==rec['rs7731706'].alleles.index('G')) for z in h])
 r=float(np.corrcoef(X,Y)[0,1]);o[pop]={'sample_total':len(people),'haplotypes_with_phased_both':len(X),'r':r,'r_squared':r*r,'freq_rs7734045_C':float(np.mean(X)),'freq_rs7731706_G':float(np.mean(Y))}
out={'source_url':url,'panel_sha256':hashlib.sha256(open('work/1kg_phase3_panel.tsv','rb').read()).hexdigest(),'variants':{k:{'position_GRCh37':p,'counted':a,'other':b,'reference_record_alleles':list(rec[k].alleles)} for k,(p,a,b) in vs.items()},'groups':o};open('outputs/stage6_hansen_chr5_reference_ld.json','w').write(json.dumps(out,indent=2)+'\n');print(json.dumps(o,indent=2))
