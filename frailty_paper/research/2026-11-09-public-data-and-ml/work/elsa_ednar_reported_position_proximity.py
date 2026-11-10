import csv,json,openpyxl
from pathlib import Path
s=openpyxl.load_workbook('work/elsa_ednar_frailty_supplement.xlsx',read_only=True,data_only=True)['Sheet1']
rows=[r for r in s.iter_rows(min_row=14,values_only=True) if isinstance(r[0],int) and isinstance(r[2],int) and isinstance(r[1],str) and r[1].startswith('rs')]
if len(rows)!=590:raise ValueError('unexpected table size')
r={}
for name in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']:
 with (Path('outputs')/name).open(newline='') as f:score=list(csv.DictReader(f))
 hits=[]
 for x in score:
  chrom=int(x.get('chr',x.get('chr_GRCh37_assumed')));pos=int(x.get('position_GRCh37_inferred',x.get('position_GRCh37_assumed')))
  for v in rows:
   d=abs(v[2]-pos)
   if chrom==v[0] and d<=1000000:hits.append({'score_variant':x['variant'],'reported_variant':v[1],'reported_distance_bp_unconfirmed_build':d})
 r[name]={'number_pairs_500kb':sum(x['reported_distance_bp_unconfirmed_build']<=500000 for x in hits),'number_pairs_1mb':len(hits),'score_variants_within_1mb':len({x['score_variant'] for x in hits}),'nearest_reported_pair':min(hits,key=lambda x:x['reported_distance_bp_unconfirmed_build']) if hits else None}
out={'source':'Mekli 2015 frailty supplement, 590 records','warning':'Coordinates compared numerically only. Source genome build not established; distance is provisional, no LD estimated, does not test the score or establish a proxy.','results':r}
Path('outputs/stage6_elsa_590_proximity_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
