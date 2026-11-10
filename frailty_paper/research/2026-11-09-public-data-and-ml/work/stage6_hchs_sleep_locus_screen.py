import openpyxl,csv,json,hashlib
from pathlib import Path
f=Path('work/hchs_sleep_health_supp.xlsx');w=openpyxl.load_workbook(f,read_only=True,data_only=True);s=w['SD15'];source=[]
for i,r in enumerate(s.iter_rows(values_only=True),1):
 if i>=5 and isinstance(r[2],str) and r[3] is not None:source.append({'sleep_score':r[1],'source_snp':r[2],'chromosome':r[3],'source_position':r[4],'gene':r[5],'effect_allele':r[6],'hchs_eaf':r[12],'hchs_beta':r[13],'hchs_se':r[14],'hchs_p':r[15]})
candidates=list(csv.DictReader(open('outputs/stage6_gsa_reference_tag_screen.csv')));ids={x['variant'] for x in candidates};matches=[r for r in source if r['source_snp'] in ids];assert len(source)==12 and len(candidates)==33
out={'source_workbook_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'sheet':'SD15','source_selected_rows':len(source),'unique_candidate_rsids':len(ids),'candidate_exact_rsid_matches':matches,'limitations':'Selected 12 loci, rsID-only; no matched row does not mean null or untested; source population phenotype differs; no genetic score validation.'}
Path('outputs/stage6_hchs_sleep_selected_locus_screen.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
