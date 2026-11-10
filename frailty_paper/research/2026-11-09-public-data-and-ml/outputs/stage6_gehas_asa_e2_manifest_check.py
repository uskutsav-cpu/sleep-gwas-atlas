import csv, hashlib, io, json, zipfile
from pyliftover import LiftOver
from pathlib import Path
root=Path(__file__).resolve().parents[1];o=root/'outputs';source=root/'work/asian_screening_e2_hg38.zip'
scores=[('15',o/'unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('22',o/'unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]
data={}
for tag,path,chrom,pos in scores:
 with path.open() as f:
  for r in csv.DictReader(f):
   k=(str(r[chrom]),int(r[pos])); entry=data.setdefault(k,{'chr':k[0],'pos37':k[1],'rsid':r['variant'],'sets':[],'name_hits':[],'coordinate_hits':[]})
   assert entry['rsid']==r['variant'];entry['sets'].append(tag)
lo=LiftOver(str(root/'work/hg19ToHg38.over.chain.gz')); unmapped=[]; lifted={}
for k,v in data.items():
 m=lo.convert_coordinate('chr'+k[0],k[1]-1)
 if len(m)!=1 or m[0][2]!='+':unmapped.append(v['rsid']);continue
 key=(m[0][0].removeprefix('chr'),m[0][1]+1)
 if key in lifted:raise ValueError('duplicate lifted coordinate')
 lifted[key]=v
byname={x['rsid']:x for x in data.values()};total=0
with zipfile.ZipFile(source) as z:
 assert len(z.namelist())==1
 f=io.TextIOWrapper(z.open(z.namelist()[0]),encoding='utf-8-sig',newline='')
 for line in f:
  if line.strip()=='[Assay]':break
 for r in csv.DictReader(f):
  if r.get('IlmnID') == '[Controls]':break
  if not r.get('IlmnID'):continue
  total+=1
  info={k:r.get(k,'') for k in ['IlmnID','Name','SNP','Chr','MapInfo','GenomeBuild','RefStrand']}
  if r['Name'] in byname:byname[r['Name']]['name_hits'].append(info)
  try:key=(str(r['Chr']),int(r['MapInfo']))
  except (ValueError,TypeError):continue
  if key in lifted:lifted[key]['coordinate_hits'].append(info)
assert total==657490, f'Unexpected assay count: {total}'
rows=[]
for v in data.values():
 rows.append({'variant':v['rsid'],'chr':v['chr'],'pos37':v['pos37'],'sets':','.join(v['sets']),'name_hits':len(v['name_hits']),'coordinate_hits':len(v['coordinate_hits']),'same_name_and_coordinate':sum(x['Name']==v['rsid'] for x in v['coordinate_hits']),'matched_probe_names':';'.join(x['Name'] for x in v['coordinate_hits']),'matched_illumina_alleles':';'.join(x['SNP'] for x in v['coordinate_hits'])})
with (o/'stage6_gehas_asa_e2_manifest_check.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
counts={tag:{'variants':sum(tag in x['sets'] for x in data.values()),'coordinate_matches':sum(bool(x['coordinate_hits']) for x in data.values() if tag in x['sets']),'same_rsid_coordinate':sum(any(q['Name']==x['rsid'] for q in x['coordinate_hits']) for x in data.values() if tag in x['sets'])} for tag in ['15','22']}
(o/'stage6_gehas_asa_e2_manifest_receipt.json').write_text(json.dumps({'manifest':'Illumina ASA-24v1-0_E2 GRCh38','source_url':'https://webdata.illumina.com/downloads/productfiles/asian-screening-array/infinium-asian-screening-array-24v1-0_E2-manifest-file-csv.zip','archive_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'unmapped_or_nonforward':unmapped,'rows_scanned':total,'counts':counts,'details':list(data.values()),'limitation':'coordinate presence is not validated allele contrast, project genotype availability, or imputation quality; exact cohort array version unknown'},indent=2)+'\n')
print(json.dumps({'unmapped_or_nonforward':unmapped,'rows_scanned':total,'counts':counts},indent=2))
