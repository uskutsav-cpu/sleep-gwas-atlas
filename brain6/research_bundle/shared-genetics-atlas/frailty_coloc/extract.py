import csv,gzip,pathlib,collections,json,hashlib,sys
kind=sys.argv[1];dest=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/frailty_coloc');dest.mkdir(exist_ok=True,parents=True)
if kind=='FI':
 path=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz'); regions={'3':(117142005,118142005),'7':(113571035,114571035)};keys=('CHR','BP')
elif kind=='R13':
 path=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz');regions={'3':(117000000,118600000),'7':(113900000,115100000)};keys=('#chrom','pos')
else:raise ValueError(kind)
rows=[];count=0
with gzip.open(path,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  count+=1;chrom=r[keys[0]]
  if chrom not in regions:continue
  pos=int(r[keys[1]]);a,b=regions[chrom]
  if a<=pos<=b:rows.append(r)
with (dest/(kind+'_extracted.tsv')).open('w',newline='') as f:
 if rows:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
print(kind,'source_rows',count,'region_rows',len(rows),'path',path,flush=True)
(dest/(kind+'_extraction.json')).write_text(json.dumps({'kind':kind,'path':str(path),'source_size':path.stat().st_size,'n_source_rows':count,'n_extracted':len(rows),'regions':regions,'extraction_sha256':hashlib.sha256((dest/(kind+'_extracted.tsv')).read_bytes()).hexdigest()},indent=2)+'\n')
