import csv,json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'target_filtered.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
a=json.loads((p/'receipt.json').read_text());assert len(r)==a['target_rows']==131
assert all(x['target']=='chr5:104484889-104484891:cg24905254' and x['study']=='Knight-ADRC' and x['context']=='PC' and x['xQTLtype']=='mQTL' for x in r)
lead=[x for x in r if x['rsID']=='rs1442112'];target=[x for x in r if x['rsID']=='rs12658032'];assert len(lead)==len(target)==1
assert float(lead[0]['pvalue'])==min(float(x['pvalue']) for x in r)
original=next(csv.DictReader((p.parent/'source.tsv').open(),delimiter='\t'));assert target[0]==original
with (p/'pairwise_ld.tsv').open() as f:ld=list(csv.DictReader(f,delimiter='\t'))
assert len(ld)==6 and len([x for x in ld if x['variant_1']=='rs1442112'])==3
assert hashlib.sha256((p/'target_filtered.tsv').read_bytes()).hexdigest()==a['output_sha256']
print('131 filtered rows; lead, exact target, and six LD pairs verified')
