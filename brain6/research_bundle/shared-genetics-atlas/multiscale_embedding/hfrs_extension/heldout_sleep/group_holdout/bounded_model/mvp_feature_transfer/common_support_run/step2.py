import csv,gzip,hashlib,json,math,os,subprocess,sys
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas'); ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');m=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged');f=Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'); feats=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];files={x:m/(x+'.sumstats.gz') for x in feats}|{x:f/(x+'.sumstats.gz') for x in targets};comp=str.maketrans('ACGT','TGCA')
def key(a,b):return min(''.join(sorted((a,b))),''.join(sorted((a.translate(comp),b.translate(comp)))))
def sha(file):
 h=hashlib.sha256()
 with file.open('rb') as fp:
  for chunk in iter(lambda:fp.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
with gzip.open(p/'nine_input_support.tsv.gz','rt') as h:common=dict(line.rstrip().split('\t') for line in h)
assert len(common)==1018491
allowed=set(); ldhash={}
for ch in range(1,23):
 path=ref/(str(ch)+'.l2.ldscore.gz');ldhash[str(ch)]=sha(path)
 with gzip.open(path,'rt',newline='') as h:
  for row in csv.DictReader(h,delimiter='\t'):
   rs=row['SNP']
   if rs in common:
    assert rs not in allowed and math.isfinite(float(row['L2']));allowed.add(rs)
common={rs:v for rs,v in common.items() if rs in allowed};print('nine inputs + LD reference',len(common),flush=True)
receipt=p/'mvp_stream_receipt.json';assert not receipt.exists();env=os.environ.copy();env['STREAM_RECEIPT']=str(receipt);source=atlas/'replication/mvp_finngen_genomewide/stream.py'; proc=subprocess.Popen(['/opt/miniconda3/bin/python3',str(source)],stdout=subprocess.PIPE,text=True,env=env,bufsize=1024*1024)
rows={}; header=proc.stdout.readline().rstrip().split('\t');assert header==['SNP','A1','A2','N','Z']
for line in proc.stdout:
 a=line.rstrip().split('\t');rs=a[0]
 if rs in common and key(a[1],a[2])==common[rs]:
  assert rs not in rows;rows[rs]='\t'.join(a)+'\n'
assert proc.wait()==0 and receipt.exists();prior=atlas/'replication/mvp_finngen_genomewide/stream_receipt.json';assert json.loads(receipt.read_text())==json.loads(prior.read_text()); common={rs:common[rs] for rs in rows}; print('final common set',len(common),flush=True);ids=sorted(common);assert len(ids)>500000
with gzip.open(p/'MVP.sumstats.gz','xt') as h:
 h.write('SNP\tA1\tA2\tN\tZ\n');h.writelines(rows[rs] for rs in ids)
del rows
receipt1=json.loads((p/'receipt.json').read_text())['input_checks'];output={'MVP':sha(p/'MVP.sumstats.gz')}
for name,path in files.items():
 assert sha(path)==receipt1[name]['sha256'];rows={}
 with gzip.open(path,'rt',newline='') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   rs=r['SNP']
   if rs in common:
    assert rs not in rows and key(r['A1'],r['A2'])==common[rs];rows[rs]='\t'.join([rs,r['A1'],r['A2'],r['N'],r['Z']])+'\n'
 assert len(rows)==len(ids)
 out=p/(name+'.sumstats.gz')
 with gzip.open(out,'xt') as h:
  h.write('SNP\tA1\tA2\tN\tZ\n');h.writelines(rows[rs] for rs in ids)
 output[name]=sha(out);print(name,len(rows),flush=True)
(p/'prepared_receipt.json').write_text(json.dumps({'common_n':len(ids),'common_ids_sha256':hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest(),'prepared_sha256':output,'ldscore_sha256':ldhash,'matching':'same sorted rsids, compatible nonpalindromic allele pairs','mvp_receipt_match':True},indent=2)+'\n')
