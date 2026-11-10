import csv,math
from pathlib import Path
p=Path('work/meis1_gtex');q={}
for line in (p/'QTD000166_MEIS1AS3_65400000_67500000.tsv').read_text().splitlines():
 x=line.split('\t');q[(int(x[2]),x[3],x[4])]=x
out=[]
with (p/'Jansen_eqtl_wide_matched.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  pos=int(r['pos38']);ref=r['ref38'];alt=r['alt38']
  if not 66350000<=pos<=66750000:continue
  if {ref,alt} in ({'A','T'},{'C','G'}) and r['rsid'] not in ('rs12621789','rs6546230'):continue
  x=q[(pos,ref,alt)]
  try:
   effect=math.log(float(r['gwas_effect']));se=float(r['gwas_se']);qb=float(x[9]);qs=float(x[10]);af=float(x[12])/float(x[13])
   if not all(map(math.isfinite,[effect,se,qb,qs,af])) or se<=0 or qs<=0:continue
  except (ValueError,ZeroDivisionError):continue
  if r['gwas_effect_allele']!=alt:effect=-effect
  out.append({'id':f'2:{pos}:{ref}:{alt}','rsid':r['rsid'],'pos38':pos,'ref':ref,'alt':alt,'qtl_alt_af':af,'qtl_beta':qb,'qtl_se':qs,'Jansen_beta':effect,'Jansen_se':se})
out.sort(key=lambda x:x['pos38']);assert len(out)==len({x['id'] for x in out})
with (p/'meis1_wide_input.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print('input rows',len(out))
