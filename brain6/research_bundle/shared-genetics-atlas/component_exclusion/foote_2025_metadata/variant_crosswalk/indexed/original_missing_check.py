import csv,json,pathlib,subprocess,sys,hashlib
p=pathlib.Path(__file__).resolve().parent;rs=sys.argv[1];assert rs in {'rs60271','rs10073929','rs40465','rs4703307','rs8085338'}
with (p.parent/'selected.tsv').open() as f:data={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
x=data[rs];target=(int(x['chrom']),int(x['pos37']));url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624046/GCST90624046.tsv'
head=subprocess.check_output(['curl','-fsSI','--max-time','30',url]).decode();hdr={line.split(':',1)[0].lower():line.split(':',1)[1].strip() for line in head.splitlines() if ':' in line};n=int(hdr['content-length']);lo=0;hi=n-1;hist=[];found=[];coordinate_rows=[]
for i in range(35):
 mid=(lo+hi)//2;start=max(0,mid-5000);end=min(n-1,mid+5000);raw=subprocess.check_output(['curl','-fsSL','--max-time','45','-r',f'{start}-{end}',url]);assert len(raw)==end-start+1
 lines=raw.split(b'\n');lines=lines[(0 if start==0 else 1):(-1 if end<n-1 else len(lines))];valid=[]
 for line in lines:
  f=line.decode().rstrip('\r').split('\t')
  if len(f)!=12 or not f[0].isdigit() or not f[1].isdigit():continue
  valid.append(((int(f[0]),int(f[1])),f))
 assert valid and all(a[0]<=b[0] for a,b in zip(valid,valid[1:]))
 hist.append({'start':start,'end':end,'first':valid[0][0],'last':valid[-1][0]})
 for coord,f in valid:
  if coord==target:
   coordinate_rows.append(f)
   if f[8]==rs and {f[2],f[3]}=={x['allele_1'],x['allele_2']}:found.append(f)
 if found:break
 if target<valid[0][0]:hi=start-1
 elif target>valid[-1][0]:lo=end+1
 else:break
 if lo>hi:break
assert len(found)<=1
r={'rsid':rs,'source_url':url,'source_etag':hdr.get('etag'),'source_length':n,'status':'EXACT_POSITIVE_RECORD' if found else 'NOT_OBSERVED_IN_SAMPLED_RANGES','exact_row':found[0] if found else None,'all_rows_at_coordinate_in_sample':coordinate_rows,'range_history':hist,'selector_sha256':hashlib.sha256((p.parent/'selected.tsv').read_bytes()).hexdigest(),'limitations':'Original file not guaranteed sorted; no find does not establish original-source absence. Partial HTTP range does not verify full source checksum.'}
(p/f'{rs}_missing_original_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(rs,r['status'],'coordinate rows',len(coordinate_rows),'ranges',len(hist),flush=True)
