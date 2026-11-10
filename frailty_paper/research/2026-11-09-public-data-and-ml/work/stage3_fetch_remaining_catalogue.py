import csv,subprocess,pathlib,concurrent.futures,json,hashlib,gzip
w=list(csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')));ds=sorted({x['dataset_r7'] for x in w if x['dataset_r7']});tasks=[]
for d in ds:
 tasks += [(d,'index',f'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{d}/{d}.all.tsv.gz.tbi',pathlib.Path('work')/(d+'.all.tsv.gz.tbi')),(d,'credible_sets',f'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/{d}/{d}.credible_sets.tsv.gz',pathlib.Path('work')/(d+'.credible_sets.tsv.gz'))]
def go(t):
 d,typ,url,path=t;tmp=path.with_suffix(path.suffix+'.download');r=subprocess.run(['curl','-fsSL','--retry','2','--max-time','180',url,'-o',str(tmp)],capture_output=True,text=True)
 if r.returncode:
  if tmp.exists():tmp.unlink()
  result={'dataset':d,'kind':typ,'url':url,'status':'download_failed','curl_exit':r.returncode,'error':r.stderr.strip()[:300]}
 else:
  if typ=='credible_sets':
   try:
    with gzip.open(tmp,'rb') as f:
     while f.read(8*1024*1024):pass
   except Exception as e:
    tmp.unlink();return {'dataset':d,'kind':typ,'url':url,'status':'invalid_gzip','error':str(e)}
  tmp.replace(path);h=hashlib.sha256(path.read_bytes()).hexdigest();result={'dataset':d,'kind':typ,'url':url,'status':'downloaded','bytes':path.stat().st_size,'local_sha256':h,'provider_checksum_verified':False}
 print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(go,tasks))
pathlib.Path('outputs/stage3_remaining_catalogue_download_receipt.json').write_text(json.dumps(results,indent=2)+'\n')
