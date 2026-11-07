#!/usr/bin/env python3
import hashlib,json,random,statistics
from pathlib import Path
p=Path(__file__).resolve().parent
rows=json.loads((p/'results.json').read_text());sets=json.loads((p/'subsets.json').read_text())
wide=json.loads((p.parent/'insomnia.json').read_text());apnea=json.loads((p.parent/'sleep_apnea.json').read_text())
ids=sorted(wide['snps']);r=random.Random(4642026)
expected=[sorted(r.sample(ids,2452)) for _ in range(20)]
assert sets==expected and len(rows)==20 and [x['replicate'] for x in rows]==list(range(1,21))
assert all(len(x)==len(set(x))==2452 and set(x)<=set(ids) for x in sets)
assert all(x['fi_failed']==(x['fi_omega']<0) for x in rows)
values=[x['fi_omega'] for x in rows];target=apnea['omega'][0][0]
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
files=[p/'PROTOCOL.md',p/'make_subsets.py',p/'subsets.json',p/'run.R',p/'results.json',p.parent/'insomnia.json',p.parent/'sleep_apnea.json']
receipt={'sha256':{str(f):sha(f) for f in files},'replicates':len(rows),'negative':sum(v<0 for v in values),
 'positive':sum(v>0 for v in values),'minimum':min(values),'maximum':max(values),'median':statistics.median(values),
 'K_range':[min(x['K'] for x in rows),max(x['K'] for x in rows)],'apnea_support_K':apnea['K'],
 'apnea_support_omega':target,'random_subsets_ge_apnea':sum(v>=target for v in values)}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
