import csv,gzip,hashlib,io,json,tarfile,collections
from pathlib import Path
p=Path(__file__).resolve().parent;base=Path('/Volumes/Brain6Work/xqtl-chr3');tar=base/'snuc_cs95.tar';man=base/'snuc_cs95.manifest.txt'
def md5(f):
 h=hashlib.md5()
 with f.open('rb') as x:
  for b in iter(lambda:x.read(2**20),b''):h.update(b)
 return h.hexdigest()
assert md5(tar)=='a278bf209cf7156cf90edabf6f7c0a4f' and md5(man)=='19d539b57a3221975255144d43ebd7fe'
manifest={r['File name']:r for r in csv.DictReader(man.open(),delimiter='\t')};rows=[];selected=[]
with tarfile.open(tar) as t:
 for member in t:
  if not(member.isfile() and 'ROSMAP_CUIMC1_2_MIT/' in member.name and any('.'+c+'.' in member.name for c in ('oli','exc')) and '_3_' in member.name and member.name.endswith('.bed.gz')):continue
  blob=t.extractfile(member).read();assert member.name in manifest and hashlib.md5(blob).hexdigest()==manifest[member.name]['md5Sum'] and len(blob)==int(manifest[member.name]['File size'])
  context='oli' if '.oli.' in member.name else 'exc';selected.append({'name':member.name,'md5':hashlib.md5(blob).hexdigest(),'size':len(blob)})
  text=gzip.decompress(blob).decode();reader=csv.DictReader(io.StringIO(text),delimiter='\t')
  for r in reader:
   if r['target_ensembl_id'] not in ('ENSG00000185565','ENSG00000239268'):continue
   rows.append({'cell':context,'target_gene':r['target_gene_symbol'],'target_id':r['target_ensembl_id'],'variant':r['variant_id'],'pip':r['PIP'],'cs95':r['cs_95'],'qtl_metadata':r['target_info'],'user_input':r['user_input']})
assert len(selected)==4
fields=['cell','target_gene','target_id','variant','pip','cs95','qtl_metadata','user_input']
with (p/'study_credible_set_rows.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
ct=collections.Counter((r['cell'],r['target_id'],r['cs95']) for r in rows)
out={'archive_md5':md5(tar),'manifest_md5':md5(man),'selected_members':selected,'target_counts':{'|'.join(k):v for k,v in ct.items()},'rows_sha256':hashlib.sha256((p/'study_credible_set_rows.tsv').read_bytes()).hexdigest(),'note':'No LINC03051 target rows in selected cell and chromosome files is not proof of no association or no complete fine-mapping elsewhere.'}
(p/'fine_mapping_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
