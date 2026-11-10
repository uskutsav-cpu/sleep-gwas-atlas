import csv,hashlib,json
from pathlib import Path
root=Path('work/phecode_x_vocab');out=[];counts={}
for version in ('1.0','1.1'):
 p=root/f'PhecodeX (version {version})/phecodeX_unrolled_ICD_CM.csv'
 with p.open(newline='') as f:rows=list(csv.DictReader(f))
 codes={c:{(r['vocabulary_id'],r.get('ICD') or r.get('icd')) for r in rows if r['phecode']==c} for c in ('SS_807','SS_807.21','NS_333.2')}
 counts[version]={'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fatigue_codes':len(codes['SS_807']),'weakness_codes':len(codes['SS_807.21']),'insomnia_codes':len(codes['NS_333.2']),'weakness_codes_in_fatigue':len(codes['SS_807']&codes['SS_807.21']),'insomnia_codes_in_fatigue':len(codes['SS_807']&codes['NS_333.2'])}
 for voc,icd in sorted(set.union(*codes.values())):
  out.append({'version':version,'vocabulary':voc,'diagnosis_code':icd,'fatigue_parent':int((voc,icd) in codes['SS_807']),'weakness_child':int((voc,icd) in codes['SS_807.21']),'insomnia':int((voc,icd) in codes['NS_333.2'])})
p=Path('outputs/stage6_aou_phecodex_mapping_overlap.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
rec={'repository':'https://github.com/PheWAS/PhecodeXVocabulary','versions':counts,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'note':'Vocabulary comparison only; does not prove the deployed AoU mapping version, case/control definition, counts or patient overlap.'}
Path('outputs/stage6_aou_phecodex_mapping_overlap_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
