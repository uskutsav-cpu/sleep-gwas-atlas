import pathlib,hashlib,gzip,csv,collections,zlib,subprocess,json,math
p=pathlib.Path(__file__).resolve().parent
ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz')
url='https://zenodo.org/records/14011550/files/GWAS_Exhaustion_maf0.01_info0.9.txt.gz?download=1'; expected='d46eb916b59abfea983d2e03d9ee7b1c'
def digest(path,alg):
 h=hashlib.new(alg)
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
refs={};duplicates=set()
with gzip.open(ref,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  k=x['SNP']
  if k in refs:duplicates.add(k)
  else:refs[k]=(x['CHR'],x['BP'],frozenset((x['A1'],x['A2'])))
for k in duplicates:refs.pop(k,None)
proc=subprocess.Popen(['curl','-fsSL','--retry','2',url],stdout=subprocess.PIPE)
md5=hashlib.md5();sha=hashlib.sha256();z=zlib.decompressobj(31);buffer=b'';size=0;header=None;counts=collections.Counter();by=collections.defaultdict(collections.Counter);seen=set();dup=set()
def row(line):
 global header
 if header is None:
  header=line.decode().split('\t');assert header==['CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','BETA','SE','CHISQ','LOG10P'],header;return
 values=line.decode().split('\t');assert len(values)==len(header);ch,pos,rs,a,b,fr,beta,se,chi,logp=values;bucket=by[ch];counts['rows']+=1;bucket['rows']+=1
 try:
  freq=float(fr);vbeta=float(beta);vse=float(se);vchi=float(chi);vp=float(logp)
  if not(all(map(math.isfinite,(freq,vbeta,vse,vchi,vp))) and 0<=freq<=1 and vse>0 and vchi>=0 and vp>=0):raise ValueError()
 except ValueError:counts['invalid_numeric']+=1;bucket['invalid_numeric']+=1
 if rs not in refs:counts['not_reference_id']+=1;bucket['not_reference_id']+=1;return
 counts['reference_id_rows']+=1;bucket['reference_id_rows']+=1
 if rs in seen:dup.add(rs)
 seen.add(rs)
 c,bp,alleles=refs[rs]
 if (ch,pos)!=(c,bp):counts['coordinate_mismatch']+=1;bucket['coordinate_mismatch']+=1;return
 counts['coordinate_match']+=1;bucket['coordinate_match']+=1
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:counts['non_snv']+=1;bucket['non_snv']+=1;return
 if frozenset((a,b))!=alleles:counts['allele_mismatch']+=1;bucket['allele_mismatch']+=1;return
 counts['allele_match']+=1;bucket['allele_match']+=1
 if frozenset((a,b)) in (frozenset('AT'),frozenset('CG')):counts['palindromic']+=1;bucket['palindromic']+=1
while chunk:=proc.stdout.read(1048576):
 size+=len(chunk);md5.update(chunk);sha.update(chunk);buffer+=z.decompress(chunk)
 lines=buffer.split(b'\n');buffer=lines.pop()
 for line in lines:
  if line:row(line)
if buffer:row(buffer)
assert proc.wait()==0 and z.eof and not z.unused_data and md5.hexdigest()==expected,(proc.returncode,z.eof,md5.hexdigest())
assert counts['reference_id_rows']==counts['coordinate_match']+counts['coordinate_mismatch']
assert counts['coordinate_match']==counts['allele_match']+counts['allele_mismatch']+counts['non_snv']
columns=['chromosome','rows','reference_id_rows','not_reference_id','coordinate_match','coordinate_mismatch','allele_match','allele_mismatch','non_snv','palindromic','invalid_numeric']
with (p/'by_chromosome.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,columns,delimiter='\t');w.writeheader()
 for ch in sorted(by,key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x)):w.writerow({'chromosome':ch,**{k:by[ch][k] for k in columns if k!='chromosome'}})
receipt={'source_url':url,'publisher_md5':expected,'observed_md5':md5.hexdigest(),'source_sha256':sha.hexdigest(),'compressed_bytes':size,'gzip_complete':z.eof,'reference_sha256':digest(ref,'sha256'),'reference_unique_ids':len(refs),'reference_duplicate_ids_excluded':len(duplicates),'source_reference_duplicate_ids':len(dup),'counts':dict(counts),'output_sha256':digest(p/'by_chromosome.tsv','sha256')}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
