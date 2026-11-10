"""Match project candidates to the selected Estonian replication table by GRCh38 position."""
import csv, json, hashlib
from pathlib import Path
import openpyxl
root=Path(__file__).resolve().parents[1]
def read(name):
 with (root/name).open() as h: return list(csv.DictReader(h))
a={r['marker']:(int(r['chr38']),int(r['pos38'])) for r in read('outputs/stage5_marker_grch37_grch38_liftover_audit.csv') if r['status']=='unique_expected_positive_strand'}
for r in read('outputs/stage6_tpmi_sleep_candidate_inventory.csv'):
 if r['ensembl_liftover_check']=='unique_positive_match':a[r['variant']]=(int(float(r['chr38'])),int(float(r['pos38'])))
for variant in ('rs113851554','rs78206187'):
 m=json.loads((root/f'work/stage6_ensembl_{variant}_mapping.json').read_text())['mappings']; assert len(m)==1 and m[0]['mapped']['strand']==1
 a[variant]=(int(m[0]['mapped']['seq_region_name']),int(m[0]['mapped']['start']))
book=root/'work/sleep_preprint_supplement.xlsx'; w=openpyxl.load_workbook(book,read_only=True,data_only=True)
rows=[r for r in list(w['ST5'].values)[2:] if len(r)>=15 and isinstance(r[1],(int,float)) and isinstance(r[2],(int,float))]
out={'workbook_sha256':hashlib.sha256(book.read_bytes()).hexdigest(),'replication_table_rows':len(rows),'sets':{}}
for f in ('unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv'):
 v=read('outputs/'+f); assert all(r['variant'] in a for r in v)
 matches=[]
 for candidate in v:
  for r in rows:
   if a[candidate['variant']]==(int(r[1]),int(r[2])):
    matches.append({'candidate':candidate['variant'],'counted':candidate['counted_allele'],'other':candidate['other_allele'],'phenotype':r[0],'reported_id':r[3],'effect_allele':r[4],'other_allele':r[5],'p':r[14]})
 out['sets'][f]={'markers':len(v),'matches':matches}
(root/'outputs/stage6_estonian_coordinate_screen.json').write_text(json.dumps(out,indent=2)+'\n')
print('table rows',len(rows),'matches by set',[len(x['matches']) for x in out['sets'].values()])
