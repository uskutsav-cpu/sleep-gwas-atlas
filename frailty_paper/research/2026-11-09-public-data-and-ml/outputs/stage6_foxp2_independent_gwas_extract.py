#!/usr/bin/env python3
import gzip,hashlib,json,csv
from pathlib import Path
root=Path(__file__).resolve().parents[1];base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw');lo=113403104;hi=115403104;records=[]
for name,path,hashname,expected in [('MVP_insomnia',base/'GCST90475826.tsv.gz','sha256','46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835'),('FinnGen_HFRS',base/'GCST90707291.tsv','md5','38d9434f50a4493def9d2c9d1cc1b4b3')]:
 h=hashlib.new(hashname)
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 assert h.hexdigest()==expected,(name,h.hexdigest());inp=gzip.open(path,'rt') if str(path).endswith('.gz') else path.open();target=root/'work'/f'stage6_foxp2_{name}_region.tsv.gz';total=kept=0
 with inp as f,gzip.open(target,'wt') as o:
  reader=csv.DictReader(f,delimiter='\t');writer=csv.DictWriter(o,fieldnames=reader.fieldnames,delimiter='\t');writer.writeheader()
  for r in reader:
   total+=1
   if r['chromosome'] in ['7','chr7']:
    try:pos=int(r['base_pair_location'])
    except (ValueError,TypeError):continue
    if lo<=pos<=hi:writer.writerow(r);kept+=1
 records.append({'source':str(path),'trait':name,'checksum':h.hexdigest(),'checksum_algorithm':hashname,'total_rows':total,'region_rows':kept,'region_grch38':f'chr7:{lo}-{hi}'})
 print(name,'rows',kept,'of',total,flush=True)
(root/'outputs/stage6_foxp2_independent_gwas_extract_receipt.json').write_text(json.dumps(records,indent=2)+'\n')
