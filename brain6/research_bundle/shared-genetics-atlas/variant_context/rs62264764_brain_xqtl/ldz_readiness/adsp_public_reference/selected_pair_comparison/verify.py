#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(q):
 with q.open() as f:return list(csv.DictReader(f,delimiter='\t'))
selected=read(p/'selected.tsv');pairs=read(p/'pair_comparison.tsv');raw=read(p/'source_pairs.tsv');stream=json.loads((p/'stream_receipt.json').read_text());summary=json.loads((p/'comparison_receipt.json').read_text())
assert len(selected)==42 and len(pairs)==42*41//2==861
matched=[x for x in pairs if x['status']=='matched'];assert len(matched)==len(raw)==50
assert all(math.isclose(abs(float(x['kgp_r2'])-float(x['adsp_r2'])),float(x['absolute_difference']),abs_tol=1e-12) for x in matched)
assert sum(float(x['absolute_difference'])>=.1 for x in matched)==summary['matched_difference_ge_0.1']==8
assert sum(float(x['absolute_difference'])>=.2 for x in matched)==summary['matched_difference_ge_0.2']==3
ind=read(p.parent/'pair.tsv');anchor=[x for x in raw if x['POS1']=='117920728' and x['POS2']=='117923158'];assert len(anchor)==len(ind)==1 and anchor[0]==ind[0]
assert stream['sha256']=='e7fbc170ffdba189d844174f1658adf9e088c529a29c9e4d407c92d24c01e4e2' and stream['matches_manifest_and_previous_full_stream']
print('42 variants, 861 pairs, 50 returned; baseline pair matches independently extracted row')
(p/'verification.json').write_text(json.dumps({'selected':len(selected),'all_pairs':len(pairs),'matched':len(matched),'independent_prior_anchor_equal':True,'source_sha256':stream['sha256'],'result_sha256':hashlib.sha256((p/'pair_comparison.tsv').read_bytes()).hexdigest(),'limitation':'Arithmetic and an independently extracted anchor do not validate all filter outputs or study-matched LD.'},indent=2)+'\n')
