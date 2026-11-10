#!/usr/bin/env python3
import csv,sys,json
models={}
for label,path in [('15','outputs/unvalidated_15_variant_score_prototype.csv'),('22','outputs/unvalidated_22_variant_score_prototype.csv')]:
 for row in csv.DictReader(open(path)):
  k=row['variant']; c=row.get('chr') or row.get('chromosome'); pos=row.get('position_GRCh37_assumed') or row.get('position_GRCh37_inferred') or row.get('position')
  if not c or not pos: raise ValueError((path,row.keys()))
  if k in models and models[k]['key']!=(str(c).replace('chr',''),int(pos)): raise ValueError(k)
  models.setdefault(k,{'key':(str(c).replace('chr',''),int(pos)),'sets':[]})['sets'].append(label)
bypos={}; byname={}
for k,v in models.items(): bypos.setdefault(v['key'],[]).append(k);byname[k.lower()]=k
reader=csv.reader(sys.stdin); header=None; hits=[]; scanned=0
for row in reader:
 if row and row[0]=='IlmnID' and 'MapInfo' in row: header=row;break
if not header: raise ValueError('header missing')
ix={x:header.index(x) for x in ['IlmnID','Name','GenomeBuild','Chr','MapInfo','SNP','RefStrand']}
for row in reader:
 if row and row[0]=='[Controls]':break
 if len(row)<len(header):continue
 scanned+=1; name=row[ix['Name']].lower(); key=(row[ix['Chr']].removeprefix('chr'),row[ix['MapInfo']]); matched=set()
 if name in byname:matched.add(byname[name])
 try:matched.update(bypos.get((key[0],int(key[1])),[]))
 except ValueError:pass
 for marker in matched:
  hits.append({'marker':marker,'model_sets':','.join(models[marker]['sets']),'match_type':('name_and_position' if name==marker.lower() and key==(models[marker]['key'][0],str(models[marker]['key'][1])) else 'name' if name==marker.lower() else 'position'),'probe_name':row[ix['Name']],'probe_id':row[ix['IlmnID']],'genome_build':row[ix['GenomeBuild']],'chromosome':key[0],'position':key[1],'snp':row[ix['SNP']],'reference_strand':row[ix['RefStrand']]})
assert scanned==4284426,scanned
with open('outputs/stage6_framingham_omni5_manifest_matches.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['marker','model_sets','match_type','probe_name','probe_id','genome_build','chromosome','position','snp','reference_strand']);w.writeheader();w.writerows(sorted(hits,key=lambda x:(x['marker'],x['probe_id'])))
print(json.dumps({'scanned':scanned,'markers':len(models),'hit_markers':len(set(x['marker'] for x in hits)),'probe_rows':len(hits)}))
