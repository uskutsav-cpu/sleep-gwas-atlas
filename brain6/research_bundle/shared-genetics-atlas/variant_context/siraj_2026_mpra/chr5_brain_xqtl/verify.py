import csv,json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent;source=p/'source.tsv'
with source.open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==1
x=r[0];assert x['rsID']=='rs12658032' and x['variantId']=='chr5:104568525:A:G' and x['xQTLtype']=='mQTL' and x['context']=='PC' and x['study']=='Knight-ADRC'
beta=float(x['beta']);se=float(x['se']);z=float(x['Z']);pv=float(x['pvalue']);q=float(x['FDR'])
assert math.isclose(beta/se,z,rel_tol=1e-6)
normal_p=math.erfc(abs(z)/math.sqrt(2))
assert pv>0 and normal_p>0
assert float(x['distance_to_target'])==83635
receipt={'source_url':'https://xqtl.niagads.org/download/rsid/rs12658032','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'returned_rows':len(r),'source_filter':'portal downloadable significant associations; no returned expression rows is not evidence of no expression association or no testing','checks':'variant allele identity, beta/SE, source-reported distance; reported p differs from standard-normal calculation and source test distribution is unverified','reported_p':pv,'two_sided_normal_p_from_z':normal_p}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
