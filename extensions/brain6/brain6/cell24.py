"""Full pair x cell family over matched independent LD blocks.

Scores are supplied from reviewed cell annotations, not inferred biological
mechanisms. Missing planned cells remain in multiplicity accounting. This is
not a substitute for raw single-cell QC, MAGMA or stratified LDSC.
"""
from __future__ import annotations
import collections
from pathlib import Path
from .stats import stratified_enrichment,bh
from .artifacts import transaction
from .io import require,read_json,read_tsv,write_json,write_tsv


def run(manifest,root,name):
    c=read_json(manifest)
    require(c.get('reviewed') is True and c.get('unit')=='independent_ld_block',
            'Review independent units, score construction and control matching')
    pairs=c['expected_pairs'];cells=c['expected_cells']
    require(pairs and len(pairs)==len(set(pairs)) and cells and len(cells)==len(set(cells)),
            'Freeze distinct nonempty pair and cell universes')
    require(c.get('exchangeability_justification'),'Matched-control exchangeability requires a justification')
    permutations=int(c.get('permutations',10000));seed=int(c.get('seed',20260908))
    require(permutations>=99,'At least 99 permutations; choose adequate FDR resolution before computation')
    raw=list(read_tsv(c['scores'],['pair_id','cell_type','locus_id','score','selected','stratum']))
    groups=collections.defaultdict(list);seen=set()
    for r in raw:
        key=r['pair_id'],r['cell_type'],r['locus_id']
        require(key[0] in pairs and key[1] in cells and key[2] and r['stratum'],'Unknown/missing cell-family unit')
        require(key not in seen,'Duplicate cell/locus score');seen.add(key)
        require(r['selected'] in {'0','1'},'Selection must be binary')
        groups[key[:2]].append(r)
    rows=[]
    for pair in pairs:
        for cell in cells:
            rs=groups[(pair,cell)];r={'pair_id':pair,'cell_type':cell,'status':'MISSING_SCORES',
                'reason':'No data for planned pair/cell','effect':None,'p':None,'n_units':len(rs),
                'permutations':permutations,'seed':seed}
            if rs:
                strata=collections.defaultdict(list)
                for x in rs:strata[x['stratum']].append(x)
                if any(not any(x['selected']=='1' for x in g) or not any(x['selected']=='0' for x in g) for g in strata.values()):
                    r.update(status='NO_EXCHANGEABLE_CONTROLS',reason='Each stratum needs selected and control loci')
                else:
                    r.update(stratified_enrichment([float(x['score']) for x in rs],[x['selected']=='1' for x in rs],
                             [x['stratum'] for x in rs],permutations=permutations,seed=seed),status='TESTED',reason='')
            rows.append(r)
    for r,q in zip(rows,bh([x['p'] for x in rows])):r['family_fdr']=q
    with transaction(root,name,stage='cell_enrichment',inputs=[manifest,c['scores']],parameters=c,
                     synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'cell_enrichment.tsv',list(rows[0]),[{k:'NA' if v is None else v for k,v in r.items()} for r in rows])
        meta['scientific_status']='PASS' if all(r['status']=='TESTED' for r in rows) else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'family_size':len(rows),
            'tested':sum(r['status']=='TESTED' for r in rows),'minimum_permutation_p':1/(permutations+1),
            'note':'Permutation resolution and matched-control assumptions must be reviewed; no causal cell mechanism is certified.'})
    return Path(root)/name
