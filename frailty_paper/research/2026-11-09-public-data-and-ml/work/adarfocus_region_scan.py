import pathlib,csv,gzip,json,math,hashlib
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort')
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
sources=[('Mak_FI',root/'fi.fastGWA','CHR','POS','SNP','P',46542374,'GRCh37 assumed'),('Mak_FP',root/'fp.fastGWA','CHR','POS','SNP','P',46542374,'GRCh37 assumed'),('MVP_insomnia',base/'GCST90475826.tsv.gz','chromosome','base_pair_location','rsid','p_value',45122459,'GRCh38'),('FinnGen_insomnia',base/'finngen_R13_F5_INSOMNIA.gz','#chrom','pos','rsids','pval',45122459,'GRCh38')]
result=[];records=[]
for name,path,chrcol,poscol,rscol,pcol,center,build in sources:
 best=[];n=0;matched=[];valid=0
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:
  hdr=f.readline().rstrip('\n\r').split('\t');ic,ip,ir,iv=[hdr.index(s) for s in (chrcol,poscol,rscol,pcol)]
  for line in f:
   n+=1; v=line.rstrip('\n\r').split('\t')
   if v[ic].removeprefix('chr')!='21':continue
   try:
    p=int(v[ip])
    if abs(p-center)>250000:continue
    prob=float(v[iv]); assert math.isfinite(prob) and 0<=prob<=1
   except (ValueError,AssertionError,IndexError):continue
   valid+=1; row={'source':name,'build':build,'chromosome':'21','position':p,'variant_id':v[ir],'p':prob,'distance_bp':p-center}
   records.append(row)
   if p==center and 'rs2838788' in v[ir].split(','):matched.append(row)
   if len(best)<5:best.append(row);best.sort(key=lambda x:x['p'])
   elif prob<best[-1]['p']:best[-1]=row;best.sort(key=lambda x:x['p'])
 result.append({'source':name,'build':build,'center':center,'window_bp':250000,'n_scanned':n,'n_valid_window_rows':valid,'candidate_hits':matched,'five_smallest_p_rows':best})
 print(name,n,valid,matched,best[:1],flush=True)
out=pathlib.Path('outputs/ad ar_region_marginal.csv'.replace(' ',''))
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
receipt={'results':result,'output_csv':out.name,'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'limitations':'Marginal p ranking, not conditional signals, fine-mapping, colocalization or common-variant comparison. Windows on different builds are candidate-centered, not identical interval liftover.'}
pathlib.Path('outputs/ad ar_region_marginal.json'.replace(' ','')).write_text(json.dumps(receipt,indent=2)+'\n')
