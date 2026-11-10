#!/usr/bin/env python3
import csv,gzip,hashlib,json
from collections import defaultdict
models={}
for label,path in [('15','outputs/unvalidated_15_variant_score_prototype.csv'),('22','outputs/unvalidated_22_variant_score_prototype.csv')]:
 for row in csv.DictReader(open(path)):
  name=row['variant']; key=(row['chr'],int(row.get('position_GRCh37_assumed') or row.get('position_GRCh37_inferred')))
  if name in models: assert models[name]['key']==key
  models.setdefault(name,{'key':key,'sets':[]})['sets'].append(label)
pos=defaultdict(list)
for name,v in models.items():pos[v['key']].append(name)
hits=[];count=0;source='work/framingham_omni5_na12878_nexus.txt.gz'
with gzip.open(source,'rt',newline='') as f:
 while f.readline().strip()!='[Data]':pass
 for row in csv.DictReader(f,delimiter='\t'):
  count+=1
  try: key=(row['Chr'],int(row['Position']))
  except ValueError:continue
  for name in pos.get(key,[]):
   hits.append({'marker':name,'score_sets':','.join(models[name]['sets']),'probe_name':row['SNP Name'],'chr':key[0],'position':key[1],'probe_snp_top_strand':row['SNP'],'name_match':str(row['SNP Name'].lower()==name.lower())})
assert count==4301332,count
with open('outputs/stage6_framingham_omni5_nexus_coordinate_matches.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['marker','score_sets','probe_name','chr','position','probe_snp_top_strand','name_match']);w.writeheader();w.writerows(sorted(hits,key=lambda a:(a['marker'],a['probe_name'])))
summary={}
for score in ('15','22'):
 selected={k for k,v in models.items() if score in v['sets']}; detected={h['marker'] for h in hits}&selected
 summary[score]={'total':len(selected),'coordinate_matches':len(detected),'matched_markers':sorted(detected),'unmatched_markers':sorted(selected-detected)}
hashobj=hashlib.sha256()
with open(source,'rb') as f:
 for block in iter(lambda:f.read(1024*1024),b''):hashobj.update(block)
summary['sha256']=hashobj.hexdigest();summary['rows']=count
with open('outputs/stage6_framingham_omni5_nexus_receipt.json','w') as f:json.dump(summary,f,indent=2);f.write('\n')
print(json.dumps(summary,indent=2))
