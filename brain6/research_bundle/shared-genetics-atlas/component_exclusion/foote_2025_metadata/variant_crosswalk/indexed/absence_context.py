import csv,json,hashlib,pathlib,sys,pysam,bisect,collections
p=pathlib.Path(__file__).resolve().parent;acc=sys.argv[1]; source=p/f'{acc}_extended.tsv';meta=p/f'{acc}.h.tsv.gz-meta.yaml'
with source.open() as f:rows=[x for x in csv.DictReader(f,delimiter='\t') if x['status']=='ABSENT']
regions=[]
for r in sorted(rows,key=lambda r:(int(r['chr38']),int(r['pos38']))):
 ch=r['chr38'];pos=int(r['pos38']);lo=max(1,pos-250000);hi=pos+250000
 if regions and regions[-1]['chrom']==ch and lo<=regions[-1]['hi']+1:regions[-1]['hi']=max(regions[-1]['hi'],hi);regions[-1]['ids'].append(r['selected_rsid'])
 else:regions.append({'chrom':ch,'lo':lo,'hi':hi,'ids':[r['selected_rsid']]})
url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/harmonised/{acc}.h.tsv.gz';tab=pysam.TabixFile(url);positions={}
for i,g in enumerate(regions):
 pos=[]
 for line in tab.fetch(g['chrom'],g['lo']-1,g['hi']):
  f=line.split('\t');assert len(f)==16 and f[0]==g['chrom'];loc=int(f[1]);assert g['lo']<=loc<=g['hi'];pos.append(loc)
 pos.sort();positions[i]=pos
byid={id:i for i,g in enumerate(regions) for id in g['ids']};out=[]
for r in rows:
 g=regions[byid[r['selected_rsid']]];pos=int(r['pos38']);a=positions[byid[r['selected_rsid']]];lo=max(1,pos-250000);hi=pos+250000;left=bisect.bisect_left(a,lo);right=bisect.bisect_right(a,hi);idx=bisect.bisect_left(a,pos);near=[abs(a[j]-pos) for j in (idx-1,idx) if 0<=j<len(a) and lo<=a[j]<=hi];out.append({'accession':acc,'rsid':r['selected_rsid'],'chr38':r['chr38'],'pos38':pos,'query_start':lo,'query_end':hi,'rows_in_250kb':right-left,'nearest_returned_bp':min(near) if near else '', 'selection':r['selection']})
with (p/f'{acc}_absence_context.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['accession','rsid','chr38','pos38','query_start','query_end','rows_in_250kb','nearest_returned_bp','selection'],delimiter='\t');w.writeheader();w.writerows(out)
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'source':url,'input_sha256':sha(source),'metadata_sha256':sha(meta),'script_sha256':sha(pathlib.Path(__file__)),'queries':[{k:v for k,v in g.items() if k!='ids'} for g in regions],'absent':len(rows),'empty_neighborhoods':sum(x['rows_in_250kb']==0 for x in out),'note':'Remote indexed response; full source checksum not verified. Cannot identify cause of gaps.'};(p/f'{acc}_absence_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(acc,'absent',len(rows),'empty ±250kb',receipt['empty_neighborhoods'],flush=True)
