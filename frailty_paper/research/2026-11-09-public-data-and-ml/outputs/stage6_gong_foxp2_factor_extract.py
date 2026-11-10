#!/usr/bin/env python3
"""Stream full Zip member, find exact GWAS IDs and allele pairs in fixed FOXP2 reference window."""
import sys,csv,json,requests,zipfile,io,hashlib
from pathlib import Path
sys.path.insert(0,'work');from stage6_remote_zip import HTTPFile
O=Path('outputs');ref=list(csv.DictReader(open(O/'stage6_foxp2_gnomad_narrow_variants.csv')));targets={r['rs']:r for r in ref};assert len(targets)==1257;j=requests.get('https://zenodo.org/api/records/19012445',timeout=30).json();f=next(a for a in j['files'] if a['key']=='MSKdg.zip');z=zipfile.ZipFile(HTTPFile(f['links']['self'],f['size']));info=z.getinfo('MSKdg.txt');found=[];n=0
with z.open(info) as raw:
 reader=csv.DictReader(io.TextIOWrapper(raw,encoding='utf-8',newline=''),delimiter='\t');assert reader.fieldnames==['SNP','CHR','BP','A1','A2','EAF','BETA','SE','Z','P','N']
 for x in reader:
  n+=1;r=targets.get(x['SNP'])
  if r and x['CHR']=='7' and x['BP']==r['pos37'] and {x['A1'],x['A2']}=={r['anchor_effect_allele'],r['anchor_other_allele']}:
   factor=1 if x['A1']==r['anchor_effect_allele'] else -1; x['rs']=x['SNP'];x['factor_beta_anchor']=factor*float(x['BETA']);x['factor_z_anchor']=factor*float(x['Z']);x['factor_eaf_anchor']=float(x['EAF']) if factor==1 else 1-float(x['EAF']);found.append(x)
  if n%1000000==0:print('scanned',n,'matches',len(found),flush=True)
assert n==4367175 and len(set(x['rs'] for x in found))==len(found)
keys=['rs','CHR','BP','A1','A2','EAF','BETA','SE','Z','P','N','factor_beta_anchor','factor_z_anchor','factor_eaf_anchor'];
with open(O/'stage6_gong_foxp2_factor_matched.csv','w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(found)
rec={'source':f['links']['self'],'published_zip_md5':f.get('checksum'),'full_zip_md5_verified':False,'full_uncompressed_member_crc_verified':True,'rows_scanned':n,'reference_input':len(ref),'matches':len(found),'source_entry_crc32':f'{info.CRC:08x}','output_sha256':hashlib.sha256((O/'stage6_gong_foxp2_factor_matched.csv').read_bytes()).hexdigest()};(O/'stage6_gong_foxp2_factor_extract_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2),flush=True)
