#!/usr/bin/env python3
"""Stream Zenodo ZIP through bounded HTTP range requests and audit all records."""
import sys,json,requests,zipfile,io,hashlib,csv
from pathlib import Path
sys.path.insert(0,'work');from stage6_remote_zip import HTTPFile
O=Path('outputs');targets=list(csv.DictReader(open(O/'unvalidated_15_variant_score_prototype.csv')));byid={r['variant']:r for r in targets};found={k:[] for k in byid};j=requests.get('https://zenodo.org/api/records/19012445',timeout=30).json();f=next(a for a in j['files'] if a['key']=='MSKdg.zip');Z=zipfile.ZipFile(HTTPFile(f['links']['self'],f['size']));info=Z.getinfo('MSKdg.txt');n=0
with Z.open(info) as raw:
 t=io.TextIOWrapper(raw,encoding='utf-8',newline='');reader=csv.DictReader(t,delimiter='\t');assert reader.fieldnames==['SNP','CHR','BP','A1','A2','EAF','BETA','SE','Z','P','N']
 for x in reader:
  n+=1
  if x['SNP'] in found:found[x['SNP']].append(x)
  if n%500000==0:print('scanned',n,flush=True)
rows=[]
for r in targets:
 candidates=found[r['variant']];matches=[x for x in candidates if x['CHR']==r['chr'] and x['BP']==r['position_GRCh37_inferred'] and {x['A1'],x['A2']}=={r['counted_allele'],r['other_allele']}];out={'variant':r['variant'],'chr37':r['chr'],'pos37':r['position_GRCh37_inferred'],'counted_allele':r['counted_allele'],'other_allele':r['other_allele'],'rsid_rows':len(candidates),'exact_rows':len(matches)}
 if len(matches)==1:
  x=matches[0];out.update(source_A1=x['A1'],source_A2=x['A2'],beta_counted=float(x['BETA'])*(1 if x['A1']==r['counted_allele'] else -1),se=float(x['SE']),p=float(x['P']),source_effect_frequency=float(x['EAF']),n_effective=x['N'])
 rows.append(out)
keys=['variant','chr37','pos37','counted_allele','other_allele','rsid_rows','exact_rows','source_A1','source_A2','beta_counted','se','p','source_effect_frequency','n_effective'];
with open(O/'stage6_gong_mskdg_15_candidate_lookup.csv','w',newline='') as h:
 writer=csv.DictWriter(h,fieldnames=keys);writer.writeheader();writer.writerows(rows)
receipt={'metadata_api':'https://zenodo.org/api/records/19012445','zip_url':f['links']['self'],'zip_bytes':f['size'],'zenodo_zip_checksum':f.get('checksum'),'whole_zip_checksum_verified':False,'zip_entry':info.filename,'zip_entry_uncompressed_size':info.file_size,'zip_entry_crc32_hex':f'{info.CRC:08x}','zip_entry_crc_verified_by_python_zipfile':True,'rows_scanned':n,'candidates':len(targets),'matched':sum(r['exact_rows']==1 for r in rows),'effective_n':'Genomic SEM working effective N, not unique participant count','limitations':'Selected candidates, prior discovery overlap with UKB and FinnGen, latent factor is not FI or FFS; pairwise source allele matching not sufficient to resolve every multiallelic representation.'};(O/'stage6_gong_mskdg_15_candidate_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
