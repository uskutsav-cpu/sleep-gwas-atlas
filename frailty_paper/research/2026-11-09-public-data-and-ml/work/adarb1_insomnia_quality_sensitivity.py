import sys,csv,gzip,math,heapq,json,pathlib,hashlib
which=sys.argv[1];base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
if which=='MVP':
 p=base/'GCST90475826.tsv.gz'; chrcol='chromosome';poscol='base_pair_location';pcol='p_value';rscol='rsid';afcol='effect_allele_frequency';qualcol='r2';alleles=('effect_allele','other_allele')
elif which=='FinnGen':
 p=base/'finngen_R13_F5_INSOMNIA.gz';chrcol='#chrom';poscol='pos';pcol='pval';rscol='rsids';afcol='af_alt';qualcol=None;alleles=('alt','ref')
else:raise ValueError(which)
center=45122459;results=[];total=0; bad=0
with gzip.open(p,'rt') as f:
 header=f.readline().rstrip('\r\n').split('\t');ix={x:header.index(x) for x in (chrcol,poscol,pcol,rscol,afcol,*alleles)}
 if qualcol:ix[qualcol]=header.index(qualcol)
 for line in f:
  v=line.rstrip('\r\n').split('\t')
  if v[ix[chrcol]].removeprefix('chr')!='21':continue
  try:
   pos=int(v[ix[poscol]])
   if abs(pos-center)>250000:continue
   total+=1;prob=float(v[ix[pcol]]);freq=float(v[ix[afcol]]);q=float(v[ix[qualcol]]) if qualcol else None
   if not (math.isfinite(prob) and 0<=prob<=1 and math.isfinite(freq) and 0<=freq<=1 and (q is None or math.isfinite(q))):raise ValueError()
  except (ValueError,IndexError):bad+=1;continue
  row={'variant':v[ix[rscol]],'position_GRCh38':pos,'allele_freq_reported':freq,'effect_allele':v[ix[alleles[0]]],'other_allele':v[ix[alleles[1]]],'p':prob,'imputation_r2':q}
  if min(freq,1-freq)>=.01 and (q is None or q>=.8):results.append(row)
results.sort(key=lambda x:x['p']);out={'source':which,'window_chr':'21','center_GRCh38':center,'radius_bp':250000,'window_rows':total,'malformed_missing_or_nonfinite_rows':bad,'passing_rows':len(results),'threshold': 'MAF >= 0.01 and MVP r2 >= 0.8' if qualcol else 'MAF >= 0.01; no imputation-quality filter in FinnGen released fields','top_ten':results[:10],'candidate_rows':[x for x in results if x['variant']=='rs2838788'],'notes':'Exploratory, post hoc unequal filter, no harmonized input support or replication; single-marker p ranking, not conditional signal or causality.'}
pout=pathlib.Path('outputs/adarb1_'+which.lower()+'_frequency_quality_sensitivity.json');pout.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
