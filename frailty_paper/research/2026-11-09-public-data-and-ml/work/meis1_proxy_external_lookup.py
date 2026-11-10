import csv,gzip,json,sys
from pathlib import Path
b=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
name=sys.argv[1];file,col={'FinnGen_insomnia':('finngen_R13_F5_INSOMNIA.gz','rsids'),'MVP_insomnia':('GCST90475826.tsv.gz','rsid'),'FinnGen_HFRS':('GCST90707291.tsv','rs_id')}[name];p=b/file;op=gzip.open if file.endswith('.gz') else open
hits=[];count=0
with op(p,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  count+=1
  if set(row[col].replace(';',',').split(',')) & {'rs13033745','rs1519102'}:hits.append(row)
out={'study':name,'source':str(p),'rows_scanned':count,'hits':hits,'scope':'exploratory post-selected SNPs; outcome and cohort overlap must be assessed; not score validation'}
Path(f'outputs/stage6_meis1_proxy_{name}_lookup.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
