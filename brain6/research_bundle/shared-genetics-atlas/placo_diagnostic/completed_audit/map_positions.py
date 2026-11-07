import csv,collections,pathlib,gzip
p=pathlib.Path(__file__).resolve().parent
with (p/'candidate_statistic_sensitivity.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
ffs={x['SNP'] for x in rows if x['pair']=='FFS'}; found=collections.defaultdict(list)
source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
with source.open() as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] in ffs: found[x['SNP']].append(x)
reference=collections.defaultdict(list)
with gzip.open('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] in ffs: reference[x['SNP']].append(x)
with (p.parent/'completed_audit_fi_only/candidate_positions.tsv').open() as f: fi={x['SNP']:x for x in csv.DictReader(f,delimiter='\t')}
for x in rows:
 if x['pair']=='FI':
  y=fi[x['SNP']]; x.update(CHR=y['CHR'],BP_GRCH37=y['BP_GRCH37'],position_status=y['position_status'])
 else:
  matches=found[x['SNP']]; assert len(matches)==1 and {matches[0]['A1'],matches[0]['A2']}=={x['A1'],x['A2']},x['SNP']
  ref=reference[x['SNP']]; assert len(ref)==1 and {ref[0]['A1'],ref[0]['A2']}=={x['A1'],x['A2']} and ref[0]['CHR']==matches[0]['CHROM'] and ref[0]['BP']==matches[0]['POS'],x['SNP']
  x.update(CHR=matches[0]['CHROM'],BP_GRCH37=matches[0]['POS'],position_status='UNIQUE_RSID_ALLELE_POSITION_MATCH_FI_REFERENCE')
 x['max_p_source']=max(float(x['sleep_p']),float(x['frailty_p']))
with (p/'candidate_positions.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for pair in ('FI','FFS'):
 xs=[x for x in rows if x['pair']==pair]
 print(pair,'chromosomes',dict(collections.Counter(x['CHR'] for x in xs)))
 for chrom in sorted({int(x['CHR']) for x in xs}):
  pos=sorted(int(x['BP_GRCH37']) for x in xs if int(x['CHR'])==chrom)
  groups=[]
  for a in pos:
   if not groups or a-groups[-1][-1]>500000: groups.append([a])
   else:groups[-1].append(a)
  print(chrom, [(g[0],g[-1],len(g)) for g in groups])
