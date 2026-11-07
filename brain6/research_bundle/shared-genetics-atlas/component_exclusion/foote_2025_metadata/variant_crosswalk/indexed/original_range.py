import csv,hashlib,json,pathlib,subprocess,sys
p=pathlib.Path(__file__).resolve().parent;rs=sys.argv[1];assert rs in {'rs62264764','rs817509','rs2431108','rs7703746'}
source=p.parent/'selected.tsv'
with source.open() as f: data={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
x=data[rs];target=(int(x['chrom']),int(x['pos37']));acc='GCST90624046'
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624046/GCST90624046.tsv'
headers=subprocess.check_output(['curl','-fsSI','--max-time','20',url]).decode();header_map={line.split(':',1)[0].lower():line.split(':',1)[1].strip() for line in headers.splitlines() if ':' in line};length=int(header_map['content-length']);assert length>1
history=[];lo=0;hi=length-1;found=[]
for iteration in range(35):
 mid=(lo+hi)//2;start=max(0,mid-5000);end=min(length-1,mid+5000)
 b=subprocess.check_output(['curl','-fsSL','--max-time','45','-r',f'{start}-{end}',url]);assert len(b)==end-start+1,(len(b),end-start+1)
 lines=b.split(b'\n');lines=lines[(0 if start==0 else 1):(-1 if end<length-1 else len(lines))];valid=[]
 for line in lines:
  f=line.decode().rstrip('\r').split('\t')
  if len(f)!=12 or not f[0].isdigit() or not f[1].isdigit():continue
  valid.append(((int(f[0]),int(f[1])),f))
 assert valid and all(a[0]<=b[0] for a,b in zip(valid,valid[1:]))
 history.append({'start':start,'end':end,'first':valid[0][0],'last':valid[-1][0]})
 for coord,f in valid:
  if coord==target and f[8]==rs and {f[2],f[3]}=={x['allele_1'],x['allele_2']}:found.append(f)
 if found:break
 if target<valid[0][0]:hi=start-1
 elif target>valid[-1][0]:lo=end+1
 else:break
 if lo>hi:break
assert len(found)==1,('could not locate unique matching row; not evidence of absence',rs,len(found))
fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rs_id','n','Q_metric','Q_metric_pvalue'];z=dict(zip(fields,found[0]));z['accession']=acc
with (p/f'{rs}_original_row.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['accession']+fields,delimiter='\t');w.writeheader();w.writerow(z)
with (p/f'{acc}_anchors.tsv').open() as f: h=next(r for r in csv.DictReader(f,delimiter='\t') if r['rsid_requested']==rs)
assert h['status']=='UNIQUE_MATCH' and {z['effect_allele'],z['other_allele']}=={h['effect_allele'],h['other_allele']}
flip=1 if z['effect_allele']==h['effect_allele'] else -1
for col in ['standard_error','p_value','Q_metric','Q_metric_pvalue']:
 assert abs(float(z[col])-float(h[col])) <= 1e-12*max(1,abs(float(z[col]))),(rs,col,z[col],h[col])
assert abs(float(z['beta'])*flip-float(h['beta']))<1e-12
r={'target':rs,'url':url,'content_length':length,'etag':header_map.get('etag'),'last_modified':header_map.get('last-modified'),'source_publisher_md5':'37f7f24f76a11a99e8840687ed847034','full_file_checksum_verified':False,'range_history':history,'source_row_sha256':hashlib.sha256(('\t'.join(found[0])+'\n').encode()).hexdigest(),'aligned_beta_sign_flip':flip,'harmonized_comparison':'beta, SE, p, Q and Q p match after effect-allele orientation; not independent source data','caution':'Original release metadata not sorted; partial range lookup verifies positive rows only.'};(p/f'{rs}_range_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(rs,'verified source row; ranges',len(history),'beta flip',flip,flush=True)
