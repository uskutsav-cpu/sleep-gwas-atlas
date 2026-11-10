import csv,json,tarfile,collections,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def read(path):
 with (root/path).open() as f:return list(csv.DictReader(f))
coords={r['marker']:(int(r['chr38']),int(r['pos38'])) for r in read('outputs/stage5_marker_grch37_grch38_liftover_audit.csv') if r['status']=='unique_expected_positive_strand'}
for r in read('outputs/stage6_tpmi_sleep_candidate_inventory.csv'):
 if r['ensembl_liftover_check']=='unique_positive_match':coords[r['variant']]=(int(float(r['chr38'])),int(float(r['pos38'])))
for rs in ('rs113851554','rs78206187'):
 v=json.loads((root/f'work/stage6_ensembl_{rs}_mapping.json').read_text())['mappings'][0]['mapped'];coords[rs]=(int(v['seq_region_name']),v['start'])
candidates={}
for file in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']:
 for r in read('outputs/'+file):
  key=r['variant'];pair=frozenset([r['counted_allele'],r['other_allele']]); candidates.setdefault(key,{'coordinate':coords[key],'pairs':[]})['pairs'].append(pair)
index=collections.defaultdict(list)
for rs,v in candidates.items(): index[v['coordinate']].append((rs,v['pairs']))
matches=[];files=0; rows=0
with tarfile.open(root/'work/varshney_eqtl_credible_sets.tar.gz','r:gz') as tar:
 for member in tar:
  if not member.isfile() or not member.name.endswith('.bed'):continue
  files+=1
  for line in tar.extractfile(member):
   cols=line.decode().strip().split('\t');rows+=1
   if len(cols)<7:continue
   try: key=(int(cols[0].removeprefix('chr')),int(cols[2]))
   except ValueError:continue
   for rs,pairs in index.get(key,[]):
    chunks=cols[3].split('_');pair=frozenset(chunks[-2:]);
    matches.append({'candidate':rs,'reported_rsid':cols[6],'variant':cols[3],'pair_match':pair in pairs,'pip':float(cols[5]),'molecular_set':member.name})
out={'source':'https://zenodo.org/records/19670721','credible_archive_md5':hashlib.md5((root/'work/varshney_eqtl_credible_sets.tar.gz').read_bytes()).hexdigest(),'files_scanned':files,'variant_rows_scanned':rows,'candidate_variants':len(candidates),'matches':matches}
(root/'outputs/stage6_varshney_muscle_candidate_screen.json').write_text(json.dumps(out,indent=2)+'\n')
print('files',files,'rows',rows,'variants',len(candidates),'matches',len(matches));print(json.dumps(matches[:25],indent=2))
