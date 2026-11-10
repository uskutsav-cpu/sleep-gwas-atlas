#!/usr/bin/env python3
import csv, hashlib, json
from pathlib import Path
from collections import defaultdict
base=Path(__file__).resolve().parent.parent
manifest=base/'work/act_human660w_v1_h.csv'
req=base/'outputs/ncds_score_variant_request.csv'
models=defaultdict(set);positions={}
with req.open(newline='') as f:
 for r in csv.DictReader(f):
  models[r['variant']].add(r['score']); p=(r['chr37'],r['pos37'])
  if r['variant'] in positions and positions[r['variant']] != p: raise ValueError('inconsistent candidate coordinate')
  positions[r['variant']]=p
bycoord=defaultdict(list)
for k,p in positions.items():bycoord[p].append(k)
h=hashlib.sha256()
with manifest.open('rb') as f:
 for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
assays=[];names=defaultdict(list);coords=defaultdict(list)
with manifest.open(newline='') as f:
 for line in f:
  if line.strip()=='[Assay]':break
 else:raise ValueError('missing assay header')
 reader=csv.DictReader(f)
 for r in reader:
  if r['IlmnID'].startswith('[Controls]'):break
  assays.append(1);name=r['Name'];p=(r['Chr'],r['MapInfo'])
  d={x:r[x] for x in ['Name','IlmnID','GenomeBuild','Chr','MapInfo','SNP','IlmnStrand','RefStrand']}
  if name in models:names[name].append(d)
  for variant in bycoord.get(p,[]):coords[variant].append(d)
if len(assays)!=655214:raise ValueError(f'Unexpected assay count {len(assays)}')
results=[]
for v in sorted(models):
 for d in coords.get(v,[]):
  results.append({'variant':v,'score_memberships':'|'.join(sorted(models[v])),'candidate_chr37':positions[v][0], 'candidate_pos37':positions[v][1],'name_exact':d['Name']==v,**d})
fields=['variant','score_memberships','candidate_chr37','candidate_pos37','name_exact','Name','IlmnID','GenomeBuild','Chr','MapInfo','SNP','IlmnStrand','RefStrand']
with (base/'outputs/stage6_act_array_coordinate_matches.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(results)
report={'manifest_url':'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/human660w-quad/human660w-quad_v1_h.csv','sha256':h.hexdigest(),'assay_count':len(assays),'candidate_count':len(models),'exact_manifest_name_count':len(names),'candidate_coordinate_count':len(coords),'name_only_variants':sorted(set(names)-set(coords)),'position_only_variants':sorted(set(coords)-set(names)),'coordinate_variants':sorted(coords),'name_variants':sorted(names),'warning':'Manufacturer v1_H GRCh37.1 manifest screen; ACT array version and post-QC actual genotypes unknown. Position and names do not establish allele, strand, presence or call quality. Not imputed coverage, validation, or outcome data.'}
(base/'outputs/stage6_act_array_manifest_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
