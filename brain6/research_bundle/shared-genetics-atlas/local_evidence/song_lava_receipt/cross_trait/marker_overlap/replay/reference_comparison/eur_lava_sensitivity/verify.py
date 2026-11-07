import csv,json,math,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent
def rows(name):
 with (p/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
allrows=rows('aligned_fi.tsv');support=rows('support.tsv');excluded=rows('excluded.tsv')
e=json.loads((p/'eur_result.json').read_text());u=json.loads((p/'matched_ukb.json').read_text())
assert len(allrows)==2741 and len(support)==u['n_snps']==e['matched_markers']==2738 and len(excluded)==3
assert [x for x in allrows if x['SNP'] in {s['SNP'] for s in support}]==support
assert {x['SNP'] for x in allrows}=={x['SNP'] for x in support}|{x['SNP'] for x in excluded}
assert len({x['SNP'] for x in allrows})==2741 and e['N_fi']==u['N_fi']==175226
assert hashlib.sha256((p/'support.tsv').read_bytes()).hexdigest()==e['support_sha256']
assert math.isclose(e['delta_sumsq']/e['K_99']-e['sigma'],e['fi_omega'],abs_tol=1e-15)
assert math.isclose(e['delta_sumsq']-e['sigma']*e['K_99'],e['fi_h2_uncapped'],abs_tol=1e-15)
assert math.isclose(u['fi_delta_sumsq']/u['K']-u['fi_sigma'],u['fi_omega'],abs_tol=1e-15)
print(json.dumps({'matched':len(support),'excluded':len(excluded),'eur_K':e['K_99'],'ukb_K':u['K'],
 'eur_fi_omega':e['fi_omega'],'ukb_fi_omega':u['fi_omega'],
 'scope':'support, input order and saved arithmetic; not independent validation of reference or statistical calibration'},indent=2))
