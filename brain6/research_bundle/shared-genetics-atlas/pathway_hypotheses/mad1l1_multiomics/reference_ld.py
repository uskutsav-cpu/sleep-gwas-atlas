import pathlib,csv,json,hashlib,collections,math,random,openpyxl,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent
wb=p/'qian_supplement.xlsx';panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');chain=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
s=openpyxl.load_workbook(wb,read_only=True,data_only=True)['Table S2a'];raw=list(s.iter_rows(values_only=True));header=raw[2];matches=[dict(zip(header,r)) for r in raw[3:] if r[0]=='cg19624444' and r[2]=='ADHD'];assert len(matches)==1
ivs=matches[0]['IVs'].split(',');assert ivs==['rs56305291','rs1637744','rs1108036'];anchor='rs28528997';variants=[anchor]+ivs
with panel.open() as f:eur={x['sample'] for x in csv.DictReader(f,delimiter='\t') if x['super_pop']=='EUR'}
assert len(eur)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);lo=LiftOver(str(chain));calls={};metadata={}
for rs in variants:
 d=json.loads((p/(rs+'_ensembl.json')).read_text());maps=[x for x in d['mappings'] if x['assembly_name']=='GRCh38' and x['seq_region_name']=='7' and x['start']==x['end']];assert len(maps)==1
 pos=maps[0]['start'];m=lo.convert_coordinate('chr7',pos-1);assert len(m)==1 and m[0][0]=='chr7' and m[0][2]=='+';bp=m[0][1]+1
 recs=[x for x in vcf.fetch('7',bp-1,bp) if x.pos==bp and x.id==rs];assert len(recs)==1 and len(recs[0].alleles)==2
 v=recs[0];ref,alt=v.alleles;data={};excluded=collections.Counter()
 for person in sorted(eur & set(vcf.header.samples)):
  sample=v.samples[person];g=sample.get('GT')
  if not g or len(g)!=2 or any(x is None or x not in (0,1) for x in g) or not sample.phased:excluded['missing_unphased_other']+=1;continue
  data[person]=tuple(int(x==1) for x in g)
 calls[rs]=data;metadata[rs]={'grch38_pos':pos,'grch37_pos':bp,'ensembl_alleles':maps[0]['allele_string'],'reference_vcf_alleles':list(v.alleles),'counted_allele':alt,'people':len(data),'excluded':dict(excluded)}
result=[]
for rs in ivs:
 shared=sorted(calls[anchor].keys() & calls[rs].keys());n=2*len(shared);pairs=[(sum(calls[anchor][i]),sum(calls[rs][i]),sum(a*b for a,b in zip(calls[anchor][i],calls[rs][i]))) for i in shared]
 def corr(xs):
  a=sum(x[0] for x in xs)/(2*len(xs));b=sum(x[1] for x in xs)/(2*len(xs));c=sum(x[2] for x in xs)/(2*len(xs));den=math.sqrt(a*(1-a)*b*(1-b));return (c-a*b)/den if den else None
 r=corr(pairs);assert r is not None;rng=random.Random(20261027);boot=[]
 for _ in range(2000):
  q=corr([pairs[rng.randrange(len(pairs))] for j in pairs]);
  if q is not None:boot.append(q*q)
 assert len(boot)==2000;boot.sort();result.append({'candidate':anchor,'candidate_counted':metadata[anchor]['counted_allele'],'instrument':rs,'instrument_counted':metadata[rs]['counted_allele'],'people':len(shared),'haplotypes':n,'r':r,'r2':r*r,'r2_boot_025':boot[49],'r2_boot_975':boot[1949]})
with (p/'instrument_ld.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
receipt={'publisher_workbook_sha256':sha(wb),'panel_sha256':sha(panel),'chain_sha256':sha(chain),'ensembl_sha256':{rs:sha(p/(rs+'_ensembl.json')) for rs in variants},'reference_url':url,'variants':metadata,'bootstrap_replicates':2000,'seed':20261027,'output_sha256':sha(p/'instrument_ld.tsv'),'caveat':'reference panel LD is not study-matched; VCF alleles do not establish published methylation effect allele or colocalization'}
(p/'instrument_ld_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))
