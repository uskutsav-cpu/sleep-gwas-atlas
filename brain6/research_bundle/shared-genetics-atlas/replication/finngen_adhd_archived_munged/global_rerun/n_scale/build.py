import gzip,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
source=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz')
target=p/'synthetic_totalN.sumstats.gz'
assert not target.exists()
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
pre=sha(source);expected=json.loads((p.parent.parent/'receipt.json').read_text())['sha256'];assert pre==expected
count=0
with gzip.open(source,'rt') as inp,gzip.open(target,'wt',compresslevel=6) as out:
 header=next(inp);assert header.strip()=='SNP\tA1\tA2\tZ\tN';out.write(header)
 for line in inp:
  values=line.rstrip('\n').split('\t');assert len(values)==5 and values[4]=='21986.309212';values[4]='495052';out.write('\t'.join(values)+'\n');count+=1
assert count==1156825 and sha(source)==pre
with gzip.open(source,'rt') as a,gzip.open(target,'rt') as b:
 assert next(a)==next(b);verified=0
 for x,y in zip(a,b):
  m=x.rstrip('\n').split('\t');n=y.rstrip('\n').split('\t');assert m[:4]==n[:4] and n[4]=='495052';verified+=1
 assert verified==count and not a.read() and not b.read()
r={'source':str(source),'source_sha256':pre,'synthetic_sha256':sha(target),'rows':count,'verified_unchanged_fields':'SNP,A1,A2,Z','old_n':'21986.309212','synthetic_n':'495052','historical_counts':[5559,489493],'effective_n_formula':4*5559*489493/(5559+489493),'note':'Synthetic alternate N is a sensitivity only, not an authenticated source or corrected estimate.'}
(p/'build_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
