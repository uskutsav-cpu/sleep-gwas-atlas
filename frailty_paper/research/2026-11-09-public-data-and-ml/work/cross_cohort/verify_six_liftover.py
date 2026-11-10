#!/usr/bin/env python3
"""Check assumed GRCh37 Lane positions against GRCh38 AoU positions, via Ensembl."""
import csv,hashlib,json,subprocess,urllib.parse
from pathlib import Path
root=Path(__file__).resolve().parents[2]
def read(name):
 with (root/'outputs'/name).open(newline='') as f:return list(csv.DictReader(f))
lane={r['variant']:r for r in read('stage6_lane2019_hunt_partners_candidate_lookup.csv')}
comp=read('stage6_lane_aou_six_variant_comparison.csv')
out=[]
for r in comp:
 source=lane[r['variant']];chrom,pos=source['chr_pos_source'].split(':');target=(r['aou_grch38_chr'],int(r['aou_grch38_pos']))
 url=f'https://grch37.rest.ensembl.org/map/human/GRCh37/{chrom}:{pos}..{pos}/GRCh38?content-type=application/json'
 raw=subprocess.check_output(['curl','-fsSL','--max-time','20',url]); m=json.loads(raw).get('mappings',[])
 dest=[x['mapped'] for x in m]
 ok=len(dest)==1 and dest[0]['seq_region_name']==target[0] and dest[0]['start']==target[1] and dest[0]['end']==target[1] and dest[0]['strand']==1
 a=set(source['source_effect_other'].split('/'));b=set(r['aou_ref_alt'].split('/'))
 out.append({'variant':r['variant'],'lane_chr_pos':source['chr_pos_source'],'aou_chr_pos':':'.join(map(str,target)),'unique_positive_coordinate_match':str(ok).lower(),'number_of_mappings':len(dest),'same_base_pair':str(a==b).lower(),'lane_effect_other':source['source_effect_other'],'aou_ref_alt':r['aou_ref_alt'],'ensembl_response_sha256':hashlib.sha256(raw).hexdigest(),'ensembl_url':url})
path=root/'outputs/stage6_lane_aou_six_liftover_audit.csv'
with path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print(f"{sum(x['unique_positive_coordinate_match']=='true' for x in out)}/{len(out)} coordinates mapped uniquely on forward strand to target; {sum(x['same_base_pair']=='true' for x in out)}/{len(out)} same base pairs")
print('output sha256',hashlib.sha256(path.read_bytes()).hexdigest())
