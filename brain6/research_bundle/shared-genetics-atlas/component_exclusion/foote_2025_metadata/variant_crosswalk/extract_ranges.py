import concurrent.futures,csv,hashlib,json,pathlib,subprocess,collections,time
p=pathlib.Path(__file__).resolve().parent;acc='GCST90624046';url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/{acc}.tsv'
with (p/'selected.tsv').open() as f:selected={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
with (p.parent/'catalog_manifest.tsv').open() as f:man={x['accession']:x for x in csv.DictReader(f,delimiter='\t')}
headers=subprocess.check_output(['curl','-fsSI','--max-time','30',url]).decode();h={line.split(':',1)[0].lower():line.split(':',1)[1].strip() for line in headers.splitlines() if ':' in line};size=int(h['content-length']);etag=h['etag'];assert not etag.startswith('W/') and size>0
step=2097152;ranges=[(i,min(size-1,i+step-1)) for i in range(0,size,step)];tag=b'\nHTTPSTATUS:'; md5=hashlib.md5();sha=hashlib.sha256();carry=b'';count=0;matches=collections.defaultdict(list);header=None;bytes_seen=0
fields_expected=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rs_id','n','Q_metric','Q_metric_pvalue']
def fetch(ab):
 a,b=ab; cmd=['curl','-fsSL','--retry','0','--connect-timeout','30','--max-time','90','--max-filesize',str(step),'-H',f'If-Match: {etag}','-r',f'{a}-{b}','-w','\nHTTPSTATUS:%{http_code}',url]
 for attempt in range(3):
  q=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE);body,sep,status=q.stdout.rpartition(tag)
  if q.returncode==0 and sep and status==b'206' and len(body)==b-a+1:return body
  if attempt==2:raise RuntimeError((a,b,q.returncode,status.decode(errors='replace'),len(body),q.stderr.decode(errors='replace')[-300:]))
  time.sleep(1+attempt)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 pending={};next_submit=0
 for i in range(len(ranges)):
  while next_submit<min(len(ranges),i+6):pending[next_submit]=ex.submit(fetch,ranges[next_submit]);next_submit+=1
  data=pending.pop(i).result();bytes_seen+=len(data);md5.update(data);sha.update(data);parts=(carry+data).split(b'\n');carry=parts.pop()
  for line in parts:
   if header is None:
    header=line.decode().rstrip('\r').split('\t');assert header==fields_expected,header;continue
   vals=line.decode().rstrip('\r').split('\t');assert len(vals)==len(header),(i,count,len(vals));count+=1
   if vals[8] in selected:matches[vals[8]].append(dict(zip(header,vals)))
  if (i+1)%20==0:print('ranges',i+1,'/',len(ranges),'bytes',bytes_seen,flush=True)
assert carry==b'' and bytes_seen==size and md5.hexdigest()==man[acc]['publisher_summary_md5'],(len(carry),bytes_seen,size,md5.hexdigest())
rows=[]
for rs,s in selected.items():
 cand=matches[rs];exact=[x for x in cand if x['chromosome'].removeprefix('chr')==s['chrom'] and x['base_pair_location']==s['pos37'] and {x['effect_allele'],x['other_allele']}=={s['allele_1'],s['allele_2']}]
 status='UNIQUE_MATCH' if len(exact)==1 else 'MULTIPLE_EXACT' if len(exact)>1 else 'ID_PRESENT_MISMATCH' if cand else 'ABSENT'
 x={'accession':acc,'rsid':rs,'selection':s['selection'],'status':status,'n_id_rows':len(cand),'n_exact':len(exact)};x.update({k:exact[0][k] if len(exact)==1 else '' for k in header});rows.append(x)
with (p/f'{acc}_crosswalk.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'url':url,'etag':etag,'last_modified':h.get('last-modified'),'bytes':size,'ranges':len(ranges),'range_bytes':step,'md5':md5.hexdigest(),'publisher_md5':man[acc]['publisher_summary_md5'],'sha256':sha.hexdigest(),'source_rows':count,'selected_sha256':hashlib.sha256((p/'selected.tsv').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'status_counts':dict(collections.Counter(x['status'] for x in rows)),'note':'Full original file MD5 verified; source is Q-filtered and not independent replication.'};(p/f'{acc}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('SUCCESS',json.dumps(receipt),flush=True)
