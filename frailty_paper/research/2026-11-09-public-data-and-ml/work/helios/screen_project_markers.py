#!/usr/bin/env python3
"""Inventory exact IDs and provisional coordinates in the public HELIOS file."""
import csv, gzip, hashlib, json
from collections import defaultdict

weights='outputs/ncds_score_variant_request.csv'
source='work/helios/novel_eqtl.txt.gz'
variants={}; positions=defaultdict(set)
with open(weights,newline='') as f:
    for row in csv.DictReader(f):
        variants[row['variant']]=True
        for build in ('37','38'):
            if row['chr'+build] and row['pos'+build]:
                positions[(row['chr'+build],int(row['pos'+build]))].add((row['variant'],build))
hits=[]; count=0; genes=set(); min_p=1.0
with gzip.open(source,'rt') as f:
    header=next(f).split()
    if header != ['chr','pos','rsid','A1','A2','gene','beta','t.stat','p.value']:
        raise ValueError('Unexpected source schema')
    for line in f:
        fields=line.split(); count+=1; genes.add(fields[5]); min_p=min(min_p,float(fields[8]))
        ids=fields[2] in variants; coords=positions.get((fields[0],int(fields[1])),set())
        if ids or coords:
            hits.append({'source':dict(zip(header,fields)), 'exact_id':ids,
                         'provisional_coordinates':[{'rsid':a,'build':b} for a,b in sorted(coords)]})
receipt={'source_sha256':hashlib.sha256(open(source,'rb').read()).hexdigest(),
         'crosswalk_sha256':hashlib.sha256(open(weights,'rb').read()).hexdigest(),
         'source_rows':count,'distinct_source_genes':len(genes),'source_min_p':min_p,
         'distinct_project_rsid':len(variants),'hits':hits,
         'limitations':'Source-selected blood cis-eQTL dataset, not whole genome or all genes; provisional coordinates do not establish build or allele identity; absence is not a null result or score validation.'}
gene_file='outputs/stage4_fuma_candidate_gene_expression.csv'
with open(gene_file,newline='') as f:
    candidate_genes={row['gene_id'].split('.')[0] for row in csv.DictReader(f) if row['gene_id']}
source_genes={gene.split('.')[0] for gene in genes}
receipt['gene_inventory']={'file':gene_file,
    'sha256':hashlib.sha256(open(gene_file,'rb').read()).hexdigest(),
    'distinct_ids':len(candidate_genes),
    'overlap':sorted(candidate_genes & source_genes),
    'scope':'Existing project FUMA expression inventory; identifiers do not establish causal genes, actual source release coverage, or lack of QTLs.'}
with open('outputs/stage6_helios_novel_eqtl_candidate_screen.json','w') as f:json.dump(receipt,f,indent=2)
print(json.dumps(receipt,indent=2))
