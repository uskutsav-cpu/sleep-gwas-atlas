import csv,gzip,json,hashlib
from pathlib import Path
specs=[('15',list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))),('22',list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv'))))]
key={};out=[]
for label,rows in specs:
 for r in rows:
  ch=r.get('chr',r.get('chr_GRCh37_assumed'));pos=int(r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed')));a=r['counted_allele'].upper();b=r['other_allele'].upper();k=(ch,pos,frozenset([a,b]));key.setdefault(k,[]).append((label,r['variant'],a,b))
f=Path('work/longevity_deelen/Results_90th_percentile.txt.gz');h=hashlib.sha256()
with f.open('rb') as inp:
 for block in iter(lambda:inp.read(1024*1024),b''):h.update(block)
with gzip.open(f,'rt') as inp:
 for row in csv.DictReader(inp,delimiter='\t'):
  a=row['EA'].upper();b=row['NEA'].upper();k=(row['Chr'],int(row['Position']),frozenset([a,b]))
  for label,rs,c,other in key.get(k,[]):
   direction=1 if a==c else -1
   out.append({'score':label,'variant':rs,'chr':row['Chr'],'position_grch37':int(row['Position']),'counted_allele':c,'other_allele':other,'eaf_counted':float(row['EAF']) if direction==1 else 1-float(row['EAF']),'beta_counted':direction*float(row['Beta']),'se':float(row['SE']),'p':float(row['P-value']),'effective_n_reported':float(row['Effective_N']),'palindromic':set([a,b]) in [set(['A','T']),set(['C','G'])]})
assert len(set((r['score'],r['variant']) for r in out))==len(out)
Path('outputs/stage6_longevity_candidate_screen.json').write_text(json.dumps({'source':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008598/Results_90th_percentile.txt.gz','sha256':h.hexdigest(),'matched':out,'requested':{label:len(rows) for label,rows in specs},'caveat':'Exact GRCh37 coordinate and unordered allele match; palindromic SNPs remain strand ambiguous. Postselected; longevity is not frailty or score validation; participant/relative overlap unresolved.'},indent=2)+'\n')
print('SHA256',h.hexdigest(),'matched by score',{label:sum(r['score']==label for r in out) for label,_ in specs})
print('nominal p<.05',[(r['score'],r['variant'],r['p']) for r in out if r['p']<.05])
