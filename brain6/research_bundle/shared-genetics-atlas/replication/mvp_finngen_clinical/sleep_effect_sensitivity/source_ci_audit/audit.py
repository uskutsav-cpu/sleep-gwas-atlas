import csv,gzip,json,math,statistics
from pathlib import Path
from statistics import NormalDist
p=Path(__file__).resolve().parent
src=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90475826.tsv.gz')
bins={'1e-8_to_1e-5':[],'1e-5_to_1e-3':[],'1e-3_to_0.05':[]};eligible=0;bad=0;discord=0;records=[]
with gzip.open(src,'rt') as h:
 r=csv.DictReader(h,delimiter='\t');headers=r.fieldnames
 for i,x in enumerate(r,1):
  if i>100000:break
  try:
   odds=float(x['odds_ratio']);lo=float(x['ci_lower']);hi=float(x['ci_upper']);pv=float(x['p_value']);b=math.log(odds)
   if not (0<lo<hi and 1e-8<=pv<=.05 and abs(b)>=.005):continue
   zci=abs(b)*2*NormalDist().inv_cdf(.975)/(math.log(hi)-math.log(lo));zs=NormalDist().inv_cdf(1-pv/2)
   ratio=zs/zci;eligible+=1
   if (lo<=1<=hi)!=(pv>=.05):discord+=1
   key='1e-8_to_1e-5' if pv<1e-5 else '1e-5_to_1e-3' if pv<1e-3 else '1e-3_to_0.05';bins[key].append(ratio)
   if len(records)<12:records.append({'rsid':x['rsid'],'odds_ratio':odds,'ci_lower':lo,'ci_upper':hi,'p_value':pv,'wald_from_ci_p':math.erfc(zci/math.sqrt(2)),'z_ratio':ratio})
  except (ValueError,TypeError,ZeroDivisionError,OverflowError):bad+=1
out={'source':str(src),'source_sha256_previously_verified':'46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835','rows_scanned':min(i,100000),'columns':headers,'eligible':eligible,'parse_errors':bad,'discordance_95ci_vs_p05':discord,'bin_stats':{k:{'n':len(v),'median_z_source_over_z_ci':statistics.median(v) if v else None,'min':min(v) if v else None,'max':max(v) if v else None} for k,v in bins.items()},'first_eligible_examples':records,'note':'First rows only, not random or independent; true source estimators and CI construction not verified.'}
(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='first_eligible_examples'},indent=2))
