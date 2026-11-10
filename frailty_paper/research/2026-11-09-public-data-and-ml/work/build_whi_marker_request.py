#!/usr/bin/env python3
"""Create a coordinate-and-allele request checklist; not a genotype or score file."""
import csv, hashlib, json
from pathlib import Path
paths={'score15':Path('outputs/unvalidated_15_variant_score_prototype.csv'),'score22_count':Path('outputs/unvalidated_22_shared_direction_score.csv'),'score22_weighted':Path('outputs/unvalidated_22_variant_score_prototype.csv')}
positions={'score15':('chr','position_GRCh37_inferred'),'score22_count':('chr_GRCh37_assumed','position_GRCh37_assumed'),'score22_weighted':('chr','position_GRCh37_assumed')}
result={}
for name,path in paths.items():
    with path.open(newline='') as f:
        for row in csv.DictReader(f):
            variant=row['variant']; chrcol,poscol=positions[name]; key=(row[chrcol],row[poscol]); d=result.setdefault(variant,{'variant':variant,'chr':key[0],'position_GRCh37_project':key[1]})
            if (d['chr'],d['position_GRCh37_project'])!=key:raise ValueError(f'Coordinate mismatch: {variant}')
            if name in d:raise ValueError(f'Duplicate within source: {variant}')
            for allele in (row['counted_allele'],row['other_allele']):
                if allele not in ('A','C','G','T'):raise ValueError(f'Non-SNV allele: {variant}')
            d[name]=row
            reference=row.get('GRCh37_ensembl_reference_alleles',row.get('reference_alleles',''))
            if reference:d.setdefault('reference_sets',set()).add(reference)
fields=['variant','chr','position_GRCh37_project','whi_snp_name_expected_if_verified','reference_alleles_project','score15_counted','score15_other','score22_count_counted','score22_count_other','score22_weighted_counted','score22_weighted_other','whi_marker_present','whi_count_allele','whi_other_allele','whi_imputation_r2','whi_allele_contrast_verified','whi_notes']
rows=[]
for d in sorted(result.values(),key=lambda v:(int(v['chr']),int(v['position_GRCh37_project']))):
    out={k:'' for k in fields};out.update(variant=d['variant'],chr=d['chr'],position_GRCh37_project=d['position_GRCh37_project'],whi_snp_name_expected_if_verified=f"{d['chr']}:{d['position_GRCh37_project']}",reference_alleles_project=';'.join(sorted(d.get('reference_sets',set()))))
    for name in paths:
        if name in d:out[name+'_counted']=d[name]['counted_allele'];out[name+'_other']=d[name]['other_allele']
    rows.append(out)
assert len(rows)==33
positions_seen=[(r['chr'],r['position_GRCh37_project']) for r in rows]
if len(set(positions_seen))!=len(positions_seen):raise ValueError('Duplicate project positions require explicit reconciliation')
out=Path('outputs/whi_score_marker_checklist.csv')
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
print(json.dumps({'n_distinct_markers':len(rows),'source_sha256':{name:hashlib.file_digest(path.open('rb'),'sha256').hexdigest() for name,path in paths.items()},'output_sha256':hashlib.file_digest(out.open('rb'),'sha256').hexdigest()},indent=2))
