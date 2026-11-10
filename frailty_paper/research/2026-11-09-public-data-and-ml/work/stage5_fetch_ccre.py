import csv,json,subprocess,concurrent.futures,pathlib,collections,urllib.parse,os
rows=list(csv.DictReader(open('outputs/stage5_score_queue.csv'))); groups=collections.defaultdict(list)
for r in rows:groups[(r['chr'],int(r['position_GRCh38'])//100000)].append(r)
p=pathlib.Path('work/stage5_ccre_raw');p.mkdir(exist_ok=True)
def fetch(item):
 (chrom,window),rs=item;start=min(int(x['position_GRCh38'])-1 for x in rs);end=max(int(x['position_GRCh38'])-1+len(x['reference_allele']) for x in rs)
 url=f'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cCREregistry;chrom={chrom};start={start};end={end}'
 target=p/f'{chrom}_{window}.json'
 if target.exists():return target.name,'cached'
 try:
  q=subprocess.run(['curl','-fsS','--max-time','45',url],capture_output=True,text=True,check=True);data=json.loads(q.stdout)
  if data.get('track')!='cCREregistry' or 'cCREregistry' not in data or data.get('error'):raise ValueError('invalid response')
  tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps({'url':url,'response':data}));os.replace(tmp,target)
  return target.name,len(data['cCREregistry'])
 except Exception as e:return target.name,'ERROR '+repr(e)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for r in ex.map(fetch,groups.items()):print(*r,flush=True)
