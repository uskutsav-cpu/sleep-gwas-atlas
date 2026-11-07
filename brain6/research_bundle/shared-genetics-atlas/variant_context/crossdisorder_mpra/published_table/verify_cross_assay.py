#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'cross_assay_published.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==4 and len({x['rsid'] for x in rows})==4
both=[r for r in rows if float(r['lee_allelic_fdr'])<.05 and float(r['siraj_allelic_fdr'])<.05];assert len(both)==1 and both[0]['rsid']=='rs62369151'
x=both[0];assert x['lee_variant_class']=='MPRA-Allelic' and x['lee_element_class']=='' and float(x['siraj_a_activity_p'])==float(x['siraj_b_activity_p'])==1
ld=json.loads((p/'rs62369151_ld_receipt.json').read_text());assert ld['n_individuals']==503 and math.isclose(ld['haplotype_r']**2,ld['r_squared'],abs_tol=1e-12)
r={'rows':len(rows),'both_source_fdr_lt_005':[x['rsid']],'ld_r2_arithmetic_checked':True,'hashes':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ('cross_assay_published.tsv','rs62369151_ld_receipt.json')},'scope':'Checks joined table and arithmetic, not source test calibration, allele direction, biological replication or study-matched LD.'}
(p/'cross_assay_published_verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
