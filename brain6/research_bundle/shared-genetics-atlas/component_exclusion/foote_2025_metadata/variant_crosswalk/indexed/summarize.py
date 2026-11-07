import csv,pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent;allrows=[]
for n in range(90624046,90624053):
 a=f'GCST{n}'
 with (p/f'{a}_anchors.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
 assert len(r)==5 and len({x['rsid_requested'] for x in r})==5
 for x in r:
  if x['status']=='UNIQUE_MATCH':
   assert 0<float(x['p_value'])<=1 and 0<float(x['Q_metric_pvalue'])<=1
   x['below_published_7factor_threshold']=str(float(x['p_value'])<7.14e-9)
  else:x['below_published_7factor_threshold']='NA'
 allrows+=r
with (p/'all_anchors.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
print('records',len(allrows),'unique',sum(x['status']=='UNIQUE_MATCH' for x in allrows),'below threshold',[(x['accession'],x['rsid_requested']) for x in allrows if x['below_published_7factor_threshold']=='True'])
