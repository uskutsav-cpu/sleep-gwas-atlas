import csv,gzip,json,sys,collections,hashlib
from pathlib import Path
from pyliftover import LiftOver
name=sys.argv[1];b=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1');files={'FI':('data/raw/frailty.txt.gz','chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','p_value'),'FFS':('frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv','CHROM','POS','A1','A2','BETA','SE','P'),'Jansen':('data/raw/insomnia.txt.gz','CHR','BP','A1','A2','OR','SE','P')}
fn,cc,pc,ac,bc,ec,sc,pv=files[name];lo=LiftOver('work/hg38ToHg19.over.chain.gz');q=collections.defaultdict(list);mapping=collections.Counter()
for line in Path('work/meis1_gtex/QTD000166_MEIS1AS3_65400000_67500000.tsv').read_text().splitlines():
 x=line.split('\t');mp=lo.convert_coordinate('chr'+x[1],int(x[2])-1)
 if len(mp)!=1 or mp[0][0]!='chr2' or mp[0][2]!='+':mapping['not_unique_plus_chr2']+=1;continue
 q[('2',str(mp[0][1]+1),x[3],x[4])].append(x);mapping['unique_plus_chr2']+=1
p=b/fn;matches=[];n=0
with (gzip.open(p,'rt') if fn.endswith('.gz') else open(p)) as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n+=1
  if row[cc]!='2':continue
  a=row[ac];o=row[bc]
  for key in [(row[cc],row[pc],o,a),(row[cc],row[pc],a,o)]:
   for x in q.get(key,[]):
    if len(q[key])!=1:continue
    matches.append({'chr38':x[1],'pos38':x[2],'ref38':x[3],'alt38':x[4],'rsid':x[-1],'pos37':row[pc],'eqtl_beta_alt':x[9],'eqtl_se':x[10],'eqtl_p':x[8],'gwas_effect_allele':a,'gwas_other_allele':o,'gwas_effect':row[ec],'gwas_se':row[sc],'gwas_p':row[pv]})
key=lambda x:(x['chr38'],x['pos38'],x['ref38'],x['alt38']);counts=collections.Counter(map(key,matches));unique=[x for x in matches if counts[key(x)]==1]
out=Path(f'work/meis1_gtex/{name}_eqtl_wide_matched.tsv');fields=['chr38','pos38','ref38','alt38','rsid','pos37','eqtl_beta_alt','eqtl_se','eqtl_p','gwas_effect_allele','gwas_other_allele','gwas_effect','gwas_se','gwas_p']
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(unique)
meta={'study':name,'input':str(p),'source_rows':n,'mapping':mapping,'matched_unique':len(unique),'duplicate_gwas_matches_excluded':len(matches)-len(unique),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'limitation':'Coordinate/pair matching only, source build/allele coding and full cis coverage require audit; no coloc'}
Path(f'outputs/stage6_meis1_{name}_eqtl_wide_receipt.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
