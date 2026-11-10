#!/usr/bin/env python3
"""Extract exact-rsID matches from Lane et al. supplementary tables 7 and 8."""
import csv, hashlib
from pathlib import Path
import openpyxl
root=Path(__file__).resolve().parents[2]
sources=[('15_variant',root/'outputs/unvalidated_15_variant_score_prototype.csv'),('22_variant_count',root/'outputs/unvalidated_22_shared_direction_score.csv')]
variants={}
for score,path in sources:
    with path.open(newline='') as fh:
        for row in csv.DictReader(fh):
            variants.setdefault(row['variant'],{})[score]=(row['counted_allele'],row['other_allele'])
wpath=root/'work/lane2019/supplement.xlsx'
w=openpyxl.load_workbook(wpath,read_only=True,data_only=True)
results={}
for tab in ('TableS7','TableS8'):
    for row in w[tab].iter_rows(min_row=5,max_col=16,values_only=True):
        rs=row[0]
        if rs not in variants: continue
        if (tab,rs) in results: raise ValueError(f'duplicate: {tab} {rs}')
        results[tab,rs]=row
out=root/'outputs/stage6_lane2019_hunt_partners_candidate_lookup.csv'
fields=['variant','score_15_counted_other','score_22_count_counted_other','table_match','chr_pos_source','nearest_gene_source','source_effect_other','hunt_or_ci','hunt_p','partners_or_ci','partners_p']
with out.open('w',newline='') as fh:
    writer=csv.DictWriter(fh,fieldnames=fields); writer.writeheader()
    for rs,scores in sorted(variants.items()):
        h=results.get(('TableS7',rs)); p=results.get(('TableS8',rs))
        if bool(h)!=bool(p): raise ValueError(f'inconsistent matching: {rs}')
        if h and (h[1],h[3])!=(p[1],p[3]): raise ValueError(f'inconsistent coding: {rs}')
        writer.writerow(dict(variant=rs,score_15_counted_other='/'.join(scores.get('15_variant',())),score_22_count_counted_other='/'.join(scores.get('22_variant_count',())),table_match='yes' if h else 'no',chr_pos_source=h[1] if h else '',nearest_gene_source=h[2] if h else '',source_effect_other=h[3] if h else '',hunt_or_ci=h[13] if h else '',hunt_p=h[14] if h else '',partners_or_ci=p[9] if p else '',partners_p=p[10] if p else ''))
print('Workbook SHA256:',hashlib.sha256(wpath.read_bytes()).hexdigest())
print('Candidates:',len(variants),'matched in both:',sum(('TableS7',v) in results for v in variants))
print('Output:',out)
