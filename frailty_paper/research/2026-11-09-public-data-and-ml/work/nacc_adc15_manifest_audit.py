import csv,json,hashlib,collections
from pathlib import Path
p=Path('work/nacc_adc15_public_manifest.csv'); rows=list(csv.DictReader(p.open()));out={}
for suffix in ['.bgen.gz','.bgen.bgi.gz','.sample']:
 v=[r for r in rows if r['file_name'].endswith(suffix)]; sizes=[int(r['file_size_bytes']) for r in v]
 out[suffix]={'file_count':len(v),'smallest_bytes':min(sizes),'largest_bytes':max(sizes),'total_bytes':sum(sizes),'requester_pays':dict(collections.Counter(r['requires_requester_pays'] for r in v))}
receipt={'manifest_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'manifest_url':'https://st1.niagads.org/portal/download-public/NG00151.v1/fm','rows':len(rows),'categories':out,'all_requester_pays':dict(collections.Counter(r['requires_requester_pays'] for r in rows)),'limitations':'Public manifest only, not actual files; listed sizes or payer settings may change, AWS, storage and compute expenses still possible; authorization and phenotype intersection remain unverified.'}
Path('outputs/stage6_nacc_adc15_manifest_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
