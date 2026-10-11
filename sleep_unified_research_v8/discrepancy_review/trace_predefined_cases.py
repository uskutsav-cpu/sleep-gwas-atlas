#!/usr/bin/env python3
"""Read-only, bounded source-direction tracing. No statistical fit or correction.
Case set fixed by assignment before extracting HGA values. Positive controls
are all Table 1 loci in de Lange 2017, independent of sleep estimates.
Private derived source extracts remain on SSD.
"""
import csv,gzip,json,math,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas')
OUT=Path(__file__).resolve().parent
HGA=Path('/Volumes/Extreme SSD/codex_sleep_gwas_atlas_v7_benchmark_20261010')
DATA=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data')
CASES=[('insomnia','ibd','3786','2031'),('insomnia','crohn','3786','2029'),('longsleep','crohn','4273','2029'),('longsleep','uc','4273','2030'),('insomnia','mdd','3786','4293')]
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as h:yield from csv.DictReader(h,delimiter='\t')
def write(p,rows):
 with p.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
v4={(x['sleep_trait'],x['outcome_trait']):x for x in read(ROOT/'sleep_unified_research_v4/tables/core_native_full_precision_rg_v4.tsv')}
meta={x['id']:x for x in read(HGA/'gwasATLAS_v20191115.txt.gz')}
keys={frozenset((s,d)) for _,_,s,d in CASES}
pairs={frozenset((x['id1'],x['id2'])):x for x in read(HGA/'gwasATLAS_v20191115_GC.txt.gz') if frozenset((x['id1'],x['id2'])) in keys}
rows=[]
for sleep,trait,sid,did in CASES:
 a=v4[sleep,trait];b=pairs[frozenset((sid,did))]
 rows.append(dict(case=f'{sleep}--{trait}',v4_rg=a['rg'],v4_se=a['se'],v4_p=a['p'],v4_ci_lower=a['ci_lower_95'],v4_ci_upper=a['ci_upper_95'],hga_sleep_id=sid,hga_outcome_id=did,hga_rg=b['rg'],hga_se=b['se'],hga_p=b['p'],hga_ci_lower=float(b['rg'])-1.96*float(b['se']),hga_ci_upper=float(b['rg'])+1.96*float(b['se']),hga_sleep_source=meta[sid]['File'],hga_outcome_source=meta[did]['File'],sign_same=(float(a['rg'])*float(b['rg'])>0),sum_rg=float(a['rg'])+float(b['rg']),difference_rg=float(a['rg'])-float(b['rg']),hga_gcov_intercept=b['gcov_int'],v4_gcov_intercept=a['gencov_intercept'],adjudication='SOURCE_DIRECTION_TRACE_PENDING' if trait!='mdd' else 'CONCORDANT_BIBLIOGRAPHIC_CONTROL_NOT_EXACT_REPLICATION'))
write(OUT/'private_predefined_case_values.tsv',rows)
write(OUT/'private_hga_selected_source_metadata.tsv',[meta[i] for i in sorted({i for c in CASES for i in c[2:]},key=int)])
r=ET.parse(OUT/'delange2017_primary.xml').getroot()
t=r.find(".//table-wrap[@id='T1']")
controls=[]
for tr in t.findall('.//tbody/tr'):
 c=[' '.join(x.itertext()).strip() for x in tr]
 controls.append(dict(snp=c[0],chr=c[1],bp=c[2],risk=c[4],nonrisk=c[5],paper_p=c[7],paper_or=c[8],disease={'IBD':'ibd','CD':'crohn','UC':'uc'}[c[10]]))
assert len(controls)==25
# Exclude strand-ambiguous and indel controls before reading local raw rows.
for c in controls:
 c['eligibility']='ELIGIBLE' if len(c['risk'])==len(c['nonrisk'])==1 and {c['risk'],c['nonrisk']} not in [{'A','T'},{'C','G'}] else 'EXCLUDED_PALINDROMIC_OR_INDEL'
results=[]
for trait in ['ibd','crohn','uc']:
 cs=[c for c in controls if c['disease']==trait]
 coords={(c['chr'],c['bp']):c for c in cs}
 raw={}
 for x in read(DATA/'raw'/f'{trait}.txt.gz'):
  coord=tuple(x['MarkerName'].split('_')[0].split(':'))
  if coord in coords:raw[coord]=x
  if len(raw)==len(coords):break
 ids={c['snp'] for c in cs}
 harmonized={x['SNP']:x for x in read(DATA/'harmonized'/f'{trait}.harmonized.tsv.gz') if x['SNP'] in ids}
 munged={x['SNP']:x for x in read(DATA/'munged'/f'{trait}.sumstats.gz') if x['SNP'] in ids}
 for c in cs:
  x=raw.get((c['chr'],c['bp']),{});h=harmonized.get(c['snp'],{});m=munged.get(c['snp'],{})
  def orient(v,a1,a2):
   if c['eligibility']!='ELIGIBLE' or not v:return None
   a1,a2=a1.upper(),a2.upper()
   if (a1,a2)==(c['risk'],c['nonrisk']):return float(v)
   if (a2,a1)==(c['risk'],c['nonrisk']):return -float(v)
   return None
  z=orient(m.get('Z',''),m.get('A1',''),m.get('A2',''));beta=orient(x.get('Effect',''),x.get('Allele1',''),x.get('Allele2',''));hb=orient(h.get('BETA',''),h.get('A1',''),h.get('A2',''))
  results.append({**c,'raw_marker':x.get('MarkerName',''),'raw_a1':x.get('Allele1',''),'raw_a2':x.get('Allele2',''),'raw_effect':x.get('Effect',''),'raw_p':x.get('P.value',''),'raw_effect_oriented_to_published_risk':beta,'harmonized_a1':h.get('A1',''),'harmonized_a2':h.get('A2',''),'harmonized_beta':h.get('BETA',''),'harmonized_beta_oriented_to_risk':hb,'munged_a1':m.get('A1',''),'munged_a2':m.get('A2',''),'munged_z':m.get('Z',''),'munged_z_oriented_to_risk':z,'raw_direction_conflicts_paper_risk':beta is not None and beta<0,'harmonized_direction_conflicts_paper_risk':hb is not None and hb<0,'munged_direction_conflicts_paper_risk':z is not None and z<0})
 print(trait, 'source controls',len(cs), 'raw matched',len(raw),'harmonized matched',len(harmonized),'munged matched',len(munged),flush=True)
write(OUT/'private_delange_paper_direction_controls.tsv',results)
summary={'case_set_fixed_before_values':True,'cases':len(CASES),'opposite_point_signs':sum(not x['sign_same'] for x in rows),'paper_table1_controls':25,'eligible_nonpalindromic_snv':sum(c['eligibility']=='ELIGIBLE' for c in controls),'raw_risk_direction_conflicts':sum(x['raw_direction_conflicts_paper_risk'] for x in results),'harmonized_risk_direction_conflicts':sum(x['harmonized_direction_conflicts_paper_risk'] for x in results),'munged_risk_direction_conflicts':sum(x['munged_direction_conflicts_paper_risk'] for x in results),'no_fit_rerun':True}
(OUT/'direction_trace_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
