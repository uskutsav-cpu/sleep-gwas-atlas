import hashlib,json,collections,openpyxl
from pathlib import Path
z=Path('work/gepitel/EWAS.zip');h=hashlib.md5()
with z.open('rb') as f:
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
assert h.hexdigest()=='d18b38b6e535d2b74c73e83cd5b22be3'
w=openpyxl.load_workbook('work/gepitel/EWAS.xlsx',read_only=True,data_only=True)
s=w['sample_sheet_clean'];header=list(next(s.iter_rows(min_row=1,max_row=1,values_only=True)));r=list(s.iter_rows(min_row=2,values_only=True));fam=[line.split() for line in open('work/gepitel/GEPITEL.fam')];ids={v for row in fam for v in row[:2]}
out={'archive_md5':h.hexdigest(),'sheets':w.sheetnames,'sample_metadata_headers':header,'sample_n':len(r),'condition_counts':dict(collections.Counter(str(x[4]) for x in r)),'first_column_is_sequential_row_numbers': [x[0] for x in r]==list(range(1,len(r)+1)),'idat_matches_genotype_id':sum(str(x[1]) in ids for x in r if x[1] is not None),'numeric_first_column_matches_genotype_ids':sum(str(x[0]) in ids for x in r if x[0] is not None),'interpretation':'Numeric row-index overlap is not a documented participant crosswalk; no linkage was performed.'}
Path('outputs/stage6_gepitel_ewas_metadata_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
