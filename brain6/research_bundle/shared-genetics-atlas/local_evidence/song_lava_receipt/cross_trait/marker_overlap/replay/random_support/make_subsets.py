#!/usr/bin/env python3
import json,random
from pathlib import Path
p=Path(__file__).resolve().parent
source=json.loads((p.parent/'insomnia.json').read_text())
ids=sorted(source['snps']);assert len(ids)==len(set(ids))==2741
r=random.Random(4642026)
samples=[sorted(r.sample(ids,2452)) for _ in range(20)]
assert all(len(s)==len(set(s))==2452 for s in samples)
(p/'subsets.json').write_text(json.dumps(samples,separators=(',',':'))+'\n')
print('Generated',len(samples),'subsets')
