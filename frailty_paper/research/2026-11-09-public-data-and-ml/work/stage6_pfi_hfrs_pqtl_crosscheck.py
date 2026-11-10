import io,zipfile,json,openpyxl
from pathlib import Path
with zipfile.ZipFile('work/pfi_supplement.zip') as z: b=z.read('ACEL-24-e70144-s001.xlsx')
s=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)['Supplementary Table 8']
panel={r[4] for r in list(s.values)[3:] if isinstance(r[5],(int,float))}
assert 'Q8TEU8' in panel
w=openpyxl.load_workbook('work/hfrs_2025_tables.xlsx',read_only=True,data_only=True)
a=[r[:7] for r in list(w['ST13'].values)[2:] if r[0]=='WFIKKN2'];assert len(a)==1 and float(a[0][5])==76.4 and float(a[0][6])==23.5
b=[r[:7] for r in list(w['ST14'].values)[2:] if r[0]=='WFIKKN2'];assert not b
out={'pfi_assay_uniprot':'Q8TEU8','assay_label':'WFKN2','hfrs_pqtl_label':'WFIKKN2','hfrs_table13_row':list(a[0]),'dementia_removed_table14_row_present':bool(b),'note':'Selected reported results; no individual-level genetics, no insomnia colocalization; omission from table 14 is not evidence of null.'}
Path('outputs/stage6_pfi_hfrs_pqtl_crosscheck.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
