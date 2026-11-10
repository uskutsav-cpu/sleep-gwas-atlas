#!/usr/bin/env python3
"""Join publicly reported candidate-level results without combining effects."""
import csv, hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parents[2]
paths={'lane':root/'outputs/stage6_lane2019_hunt_partners_candidate_lookup.csv','aou':root/'outputs/stage6_aou_eur_candidate_association.csv'}
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
lane=[r for r in read(paths['lane']) if r['table_match']=='yes']
aou=read(paths['aou']); out=[]
for r in lane:
    sources={}
    for pheno in ['NS_333.2','SS_807.21']:
        hit=[x for x in aou if x['candidate']==r['variant'] and x['phenotype']==pheno]
        if len(hit)!=1 or hit[0]['status']!='matched': raise ValueError((r['variant'],pheno,len(hit)))
        sources[pheno]=hit[0]
    a,b=sources.values()
    if (a['chrom'],a['position'],a['ref'],a['alt'])!=(b['chrom'],b['position'],b['ref'],b['alt']): raise ValueError('inconsistent AoU variant representation')
    out.append({'variant':r['variant'],'lane_effect_other':r['source_effect_other'],'hunt_or_ci':r['hunt_or_ci'],'hunt_p':r['hunt_p'],'partners_or_ci':r['partners_or_ci'],'partners_p':r['partners_p'],'aou_grch38_chr':a['chrom'],'aou_grch38_pos':a['position'],'aou_ref_alt':a['ref']+'/'+a['alt'],'aou_insomnia_p':a['pvalue'],'aou_insomnia_returned_count':a['returned_count'],'aou_weakness_p':b['pvalue'],'aou_weakness_returned_count':b['returned_count']})
outpath=root/'outputs/stage6_lane_aou_six_variant_comparison.csv'
with outpath.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=out[0].keys());writer.writeheader();writer.writerows(out)
receipt={'inputs':{k:{'file':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},'output_sha256':hashlib.sha256(outpath.read_bytes()).hexdigest(),'variants':len(out),'note':'P-values only; effect estimates are not pooled or harmonized.'}
(root/'outputs/stage6_lane_aou_six_variant_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
