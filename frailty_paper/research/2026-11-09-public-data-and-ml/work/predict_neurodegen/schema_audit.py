import glob,openpyxl,hashlib,json,pathlib
result={'source':'https://zenodo.org/records/20734280','archive_md5':hashlib.md5(pathlib.Path('work/predict_neurodegen/data.zip').read_bytes()).hexdigest(),'sheets':{},'warning':'Counts nonblank, not proof of valid measurements or linkage approval; no individual rows exported.'}
assert result['archive_md5']=='c324b026822df5630c1e126cb141de17'
for path in glob.glob('work/predict_neurodegen/*.xlsx'):
 for s in openpyxl.load_workbook(path,read_only=True,data_only=True):
  it=s.iter_rows(values_only=True);headers=[str(x) if x is not None else '' for x in next(it)];counts=[0]*len(headers);ids=[];num=0
  for row in it:
   if not any(v is not None for v in row):continue
   num+=1
   for i,v in enumerate(row[:len(headers)]):
    if v is not None:counts[i]+=1
   if row[0] is not None:ids.append(str(row[0]).strip())
  result['sheets'][s.title]={'nonempty_rows':num,'rows_with_nonempty_id':len(ids),'unique_ids':len(set(ids)),'column_headers':headers,'nonblank_counts':dict(zip(headers,counts))}
pathlib.Path('outputs/stage6_predict_neurodegen_schema_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:(v['nonempty_rows'],v['unique_ids']) for k,v in result['sheets'].items()})
