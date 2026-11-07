#!/usr/bin/env python3
import hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
restricted=json.loads((p/'restricted.json').read_text())
insomnia=json.loads((p.parent/'insomnia.json').read_text())
apnea=json.loads((p.parent/'sleep_apnea.json').read_text())
assert set(restricted['snps'])==set(apnea['snps']) and len(restricted['snps'])==2452
assert restricted['K']==apnea['K']==185 and insomnia['K']==216
assert restricted['N'][0]==apnea['N'][0]==insomnia['N'][0]
for key in ('omega','sigma'):
 assert math.isclose(restricted[key][0][0],apnea[key][0][0],rel_tol=1e-12,abs_tol=1e-17)
assert math.isclose(restricted['delta_squared_sums'][0],apnea['delta_squared_sums'][0],rel_tol=1e-12)
assert math.isclose(restricted['h2_obs_capped'][0],apnea['h2_obs_capped'][0],rel_tol=1e-12)
assert restricted['failed'][0] is False and insomnia['failed'][0] is True
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
inputs=[p/'PROTOCOL.md',p/'run.R',p/'restricted.json',p.parent/'insomnia.json',p.parent/'sleep_apnea.json',p.parent/'insomnia_quoted_input.tsv']
receipt={'input_sha256':{str(f):sha(f) for f in inputs},'exact_snp_set_match':True,'same_fi_N':True,
 'wide_insomnia_fi_omega':insomnia['omega'][0][0], 'restricted_insomnia_fi_omega':restricted['omega'][0][0],
 'apnea_fi_omega':apnea['omega'][0][0], 'restricted_K':restricted['K'], 'diagnostic_only':True}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
