import csv,hashlib,json,re,xml.etree.ElementTree as ET
from pathlib import Path
p=Path('work/kripke_paper.xml');root=ET.parse(p).getroot(); ids=set(re.findall(r'\brs\d+\b',' '.join(root.itertext())))
res={'source':'https://doi.org/10.1016/j.sleep.2014.11.003','xml_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'article_distinct_rsids':len(ids),'score_intersections':{},'warning':'Only identifiers occurring in full article XML, not the complete genotyping chip or unavailable supplementary workbook. No genotype or score association was evaluated.'}
for name in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']:
 with (Path('outputs')/name).open(newline='') as f:score={r['variant'] for r in csv.DictReader(f)}
 res['score_intersections'][name]=len(ids&score)
Path('outputs/stage6_scripps_sleep_article_rsid_screen.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
