import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;atlas=p.parents[2];r=json.loads((p/'receipt.json').read_text());f=json.loads((p/'ffs_receipt.json').read_text());v=json.loads((p/'ensembl.json').read_text());prov=json.loads((atlas.parent/'brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256'];assert r['mapping38']==104568698
assert any(x.get('assembly_name')=='GRCh38' and x.get('seq_region_name')=='5' and x.get('start')==104568698 and x.get('allele_string')=='G/A/C' for x in v['mappings'])
for k in ('brain_insomnia','brain_adhd'):
 s=r['source_reports'][k];assert s['sha256']==prov[pathlib.Path(s['path']).name] and not s['rows']
assert not r['source_reports']['FI']['rows'];rr=r['source_reports']['R13']['rows'];assert len(rr)==1 and (rr[0]['ref'],rr[0]['alt'])==('G','A')
assert len(f['matches'])==1 and (f['matches'][0]['A1'],f['matches'][0]['A2'])==('G','A')
for row in csv.DictReader((p/'ld.tsv').open(),delimiter='\t'):
 a,b,c,d=(int(row[k]) for k in ('n00','n01','n10','n11'));assert a+b+c+d==2*int(row['n']);x=(a*d-b*c)**2/((a+b)*(c+d)*(a+c)*(b+d));assert math.isclose(x,float(row['r2']),abs_tol=1e-12)
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();(p/'verification.json').write_text(json.dumps(dict(status='PASS',checks='two frozen input hashes, source matches, Ensembl alleles, LD count arithmetic; bootstrap and source calibration not independently checked',ld_sha256=sha(p/'ld.tsv'),receipt_sha256=sha(p/'receipt.json'),ffs_receipt_sha256=sha(p/'ffs_receipt.json'),ensembl_sha256=sha(p/'ensembl.json')),indent=2)+'\n');print('PASS: input hashes, alleles, source rows and reference LD arithmetic')
