import csv,pysam,numpy as np,json,pathlib,hashlib
from pyliftover import LiftOver
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
sources=['outputs/adarb1_mvp_frequency_quality_sensitivity.json','outputs/adarb1_finngen_frequency_quality_sensitivity.json']
raw=[json.load(open(p))['top_ten'][0] for p in sources];rows=[{'variant':'rs2838788','position_GRCh38':45122459,'effect_allele':'T','other_allele':'C'}]+raw
chain=pathlib.Path('work/hg38ToHg19.over.chain.gz');lo=LiftOver(str(chain));pan=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in pan if r['super_pop']=='EUR'};f=pysam.VariantFile(url);data={};vrec=[]
for row in rows:
 rs=row['variant'];pos=row['position_GRCh38'];a=row['effect_allele'];b=row['other_allele'];maps=lo.convert_coordinate('chr21',pos-1) or [];plus=[m for m in maps if m[0]=='chr21' and m[2]=='+'];entry={'variant':rs,'pos38':pos,'allele':a,'other':b,'n_map':len(maps)}
 if len(maps)!=1 or len(plus)!=1:entry['status']='mapping not unique positive chr21';vrec.append(entry);continue
 p=plus[0][1]+1;entry['pos37']=p;hits=[v for v in f.fetch('21',p-1,p) if v.id==rs and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={a,b}];entry['n_matching_records']=len(hits)
 if len(hits)!=1:entry['status']='no unique exact biallelic match';vrec.append(entry);continue
 var=hits[0];ix=var.alleles.index(a);calls={}
 for s in eur & set(f.header.samples):
  gt=var.samples[s].get('GT')
  if gt and len(gt)==2 and all(z is not None for z in gt) and var.samples[s].phased:calls[s]=[int(z==ix) for z in gt]
 data[rs]=calls;flat=[z for v in calls.values() for z in v];entry.update(status='matched',n_called=len(calls),frequency=float(np.mean(flat)))
 vrec.append(entry)
pairs=[]
for row in raw:
 rs=row['variant'];res={'variant':rs}
 if rs not in data or 'rs2838788' not in data:res['status']='not calculated'
 else:
  ids=sorted(data[rs].keys() & data['rs2838788'].keys());a=np.array([z for s in ids for z in data[rs][s]]);b=np.array([z for s in ids for z in data['rs2838788'][s]])
  if not len(a) or a.std()==0 or b.std()==0:res['status']='insufficient variation'
  else:
   r=float(np.corrcoef(a,b)[0,1]);res.update(status='calculated',r=r,r2=r*r,n_individuals=len(ids),n_haplotypes=len(a))
 pairs.append(res)
out={'reference':url,'population':'503 European-labeled 1000 Genomes samples; phased haplotypes','chain_sha256':hashlib.sha256(chain.read_bytes()).hexdigest(),'source_quality_receipt_hashes':[hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() for p in sources],'variants':vrec,'pairs':pairs,'limitations':'post hoc peaks, not conditional or colocalization analysis; EUR reference is not Finnish cohort; no individual data saved'}
pathlib.Path('outputs/adarb1_filtered_peak_reference_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
