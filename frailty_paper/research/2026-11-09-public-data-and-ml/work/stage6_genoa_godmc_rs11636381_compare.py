import csv,json,hashlib
from pathlib import Path
p=Path('work/godmc_rs11636381.json');j=json.loads(p.read_text())['assoc_meta'];g=list(csv.DictReader(open('outputs/stage6_genoa_rs11636381_mqtl.csv')))
assert len(j)==29 and all(x['snp']=='chr15:74134883:SNP' and {x['a1'],x['a2']}=={'C','T'} for x in j)
assert len(g)==38 and all(x['ps']=='74134883' and {x['allele1'],x['allele0']}=={'C','T'} for x in g)
a={x['cpg'] for x in j};b={x['CpG'] for x in g};assert len(a)==29 and len(b)==38
receipt={'GoDMC_rows':len(j),'GENOA_rows':len(g),'shared_cpg':sorted(a&b),'godmc_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source':'http://api.godmc.org.uk/v0.1/assoc_meta/rsid/rs11636381','warning':'GoDMC preselected, different assay and cis windows; HTTP unverified, no overlap is not negative replication.'}
Path('outputs/stage6_genoa_godmc_rs11636381_compare.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
