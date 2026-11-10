from pathlib import Path
import csv,gzip,hashlib,json
O=Path('outputs');W=Path('work');p=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90475826.tsv.gz');h=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
assert h.hexdigest()=='46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835';out=W/'stage6_rbm6_mvp_original_region.tsv.gz';total=kept=0
with gzip.open(p,'rt') as f,gzip.open(out,'wt') as o:
 rd=csv.DictReader(f,delimiter='\t');wr=csv.DictWriter(o,fieldnames=rd.fieldnames,delimiter='\t');wr.writeheader()
 for r in rd:
  total+=1
  if r['chromosome']=='3' and r['base_pair_location'].isdigit() and 49714152<=int(r['base_pair_location'])<=50711834:wr.writerow(r);kept+=1
receipt={'source_path':str(p),'source_sha256':h.hexdigest(),'matches_previously_recorded_expected_sha256':True,'region':'GRCh38 chr3:49714152-50711834','total_rows':total,'kept':kept,'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()};(O/'stage6_rbm6_mvp_original_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
