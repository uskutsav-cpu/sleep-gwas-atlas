import csv,gzip,json,sys
from pathlib import Path
base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1')
sources={'Atkins_FI':base/'data/raw/frailty.txt.gz','Jansen_insomnia':base/'data/raw/insomnia.txt.gz','Ye_FFS':base/'frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'}
name=sys.argv[1];p=sources[name];op=gzip.open if p.suffix=='.gz' else open
variants={'rs13033745','rs1519102'}; hits=[];count=0
with op(p,'rt') as f:
 reader=csv.DictReader(f,delimiter='\t')
 for row in reader:
  count+=1
  if row.get('variant_id',row.get('SNP','')) in variants:hits.append(row)
out={'study':name,'source':str(p),'rows_scanned':count,'variants':sorted(variants),'hits':hits,'scope':'post-selection single variant lookups, not independent score validation'}
Path(f'outputs/stage6_meis1_proxy_{name}_lookup.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
