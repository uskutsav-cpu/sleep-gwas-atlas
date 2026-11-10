import csv,gzip,json,sys,collections,hashlib
from pathlib import Path
name=sys.argv[1];base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
files={'MVP':('GCST90475826.tsv.gz','rsid'),'HFRS':('GCST90707291.tsv','rs_id')};fn,rs=files[name];p=base/fn
q=collections.defaultdict(list)
for line in Path('work/meis1_gtex/QTD000166_MEIS1AS3_65400000_67500000.tsv').read_text().splitlines():
 x=line.split('\t');q[(x[1],x[2],x[3],x[4])].append(x)
seen=collections.Counter();matched=[];n=0
with (gzip.open(p,'rt') if fn.endswith('.gz') else open(p)) as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n+=1
  if row['chromosome']!='2':continue
  pos=row['base_pair_location'];a=row['effect_allele'];b=row['other_allele']
  for key,order in [((row['chromosome'],pos,b,a),'same'),((row['chromosome'],pos,a,b),'swapped')]:
   if key in q:
    for x in q[key]:
     if len(q[key])!=1:continue
     seen[key]+=1;matched.append({'chr':key[0],'pos':key[1],'ref':key[2],'alt':key[3],'eqtl_rsid':x[-1],'eqtl_beta':x[9],'eqtl_se':x[10],'eqtl_p':x[8],'gwas_effect_allele':a,'gwas_other_allele':b,'gwas_beta_or':row.get('beta',row.get('odds_ratio')),'gwas_se':row.get('standard_error'),'gwas_p':row['p_value'],'alignment':order})
counts=collections.Counter((x['chr'],x['pos'],x['ref'],x['alt']) for x in matched); unique=[x for x in matched if counts[(x['chr'],x['pos'],x['ref'],x['alt'])]==1]
out=Path(f'work/meis1_gtex/{name}_eqtl_wide_matched.tsv')
with out.open('w',newline='') as f:
 fields=['chr','pos','ref','alt','eqtl_rsid','eqtl_beta','eqtl_se','eqtl_p','gwas_effect_allele','gwas_other_allele','gwas_beta_or','gwas_se','gwas_p','alignment'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(unique)
meta={'source':str(p),'gwas_rows_scanned':n,'eqtl_rows':sum(map(len,q.values())),'exact_allele_match_unique':len(unique),'duplicate_gwas_matches_excluded':len(matched)-len(unique),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'note':'GRCh38 coordinate and allele matching only; SNP rsID duplicates and multiallelic contrasts require separate review; no colocalization performed'}
Path(f'outputs/stage6_meis1_{name}_eqtl_wide_receipt.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
