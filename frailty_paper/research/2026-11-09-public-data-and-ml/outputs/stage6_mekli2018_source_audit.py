from pathlib import Path
import openpyxl,csv,json,hashlib
O=Path('outputs');p=Path('work/stage6_mekli2018_s3.xlsx');w=openpyxl.load_workbook(p,read_only=True,data_only=True);s=w['GenAssoc'];rows=list(s.iter_rows(values_only=True));header=rows[0];data=[dict(zip(header,r)) for r in rows[1:] if r[16] is not None];assert len(data)==32
for fn in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_variant_score_prototype.csv']:
 c={r['variant'] for r in csv.DictReader((O/fn).open())};over={str(r['snp']) for r in data}&c;print(fn,'exact rsID overlaps',len(over));assert not over
lead=[r for r in data if r['snp']=='rs7134291'];assert len(lead)==1
with (O/'stage6_mekli2018_published_fi_leads.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(header));writer.writeheader();writer.writerows(data)
(O/'stage6_mekli2018_source_audit.json').write_text(json.dumps({'supplement_url':'https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0207824.s003&type=supplementary','supplement_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'selected_fi_source_rows':len(data),'exact_rsid_overlaps_15':0,'exact_rsid_overlaps_22':0,'scope':'32 selected HRS hits only; non-overlap does not imply absent in HRS or null associations'},indent=2)+'\n');print('lead',lead[0])
