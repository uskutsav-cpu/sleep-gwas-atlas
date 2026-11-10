import csv,hashlib,json,openpyxl
from pathlib import Path
p=Path('outputs'); f=Path('work/charge_sleep_supp/41467_2019_8917_MOESM9_ESM.xlsx')
ids={r['variant'] for r in csv.DictReader((p/'stage6_charge_gait_song124_lookup.csv').open())};assert len(ids)==124
s=openpyxl.load_workbook(f,read_only=True,data_only=True).active; main=[]; proxy=[]
for r in s.iter_rows(min_row=5,values_only=True):
 if isinstance(r[0],str) and r[0].startswith('rs'):main.append(r[0])
 if isinstance(r[4],str) and r[4].startswith('rs'):proxy.append(r[4])
o={'candidate_count':len(ids),'ukb_duration_lead_rows':len(main),'charge_reported_proxy_rows':len(proxy),'exact_candidate_lead_matches':sorted(ids.intersection(main)),'exact_candidate_proxy_matches':sorted(ids.intersection(proxy)),'source_workbook_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'candidate_file_sha256':hashlib.sha256((p/'stage6_charge_gait_song124_lookup.csv').read_bytes()).hexdigest(),'scope':'Selected UK Biobank sleep duration loci and CHARGE proxies only; not a complete scan, score validation, test of insomnia or full frailty.'}
(p/'stage6_charge_sleep_supp_intersection.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
