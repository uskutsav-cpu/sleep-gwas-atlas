import csv, json, math
from pathlib import Path
p=Path(__file__).resolve().parent
selected=list(csv.DictReader((p/'selected_profiles.tsv').open(),delimiter='\t'))
receipt=json.loads((p.parent/'adhd_unfiltered/receipt.json').read_text())
comparisons=[]
for row in selected:
    tag=f"fusion.{row['tissue']}.{row['model']}.PGC_ADHD_EUR_2017.tsv"
    matches=[]
    for entry in receipt['matches']:
        if entry['member'].endswith(tag):
            record=dict(zip(entry['header'].split('\t'),entry['row'].split('\t')))
            if record['ID']==row['feature']:matches.append(record)
    assert len(matches)==1
    x=matches[0]
    assert x['MODEL']==row['selected_method'] and x['NSNP']==row['reported_nsnps']
    # Some association rows round the profile R² to two decimal places.
    assert abs(float(x['MODELCV.R2'])-float(row['reported_method_r2'])) <= 0.00500001
    assert math.isclose(float(x['MODELCV.PV']),float(row['reported_method_p']),rel_tol=1e-12)
    comparisons.append({'tissue':row['tissue'],'feature':row['feature'],'model':x['MODEL'],'cv_r2':float(x['MODELCV.R2'])})
assert len(comparisons)==8
print(json.dumps({'matched':len(comparisons),'results':comparisons},indent=2))
