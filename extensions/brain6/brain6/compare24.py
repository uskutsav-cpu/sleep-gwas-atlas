"""Cross-disorder independent-block overlap with explicit coverage and matched nulls."""
from __future__ import annotations
import collections
import itertools
from pathlib import Path
from .pathway24 import conditional_overlap
from .io import require, read_tsv, read_json, write_tsv, write_json
from .stats import bh
from .artifacts import transaction


def compare(manifest, root, name):
    c=read_json(manifest)
    require(c.get('reviewed') is True and c.get('unit')=='independent_ld_block',
            'Common LD-block universe and matching must be reviewed')
    strategy=c.get('dependence_strategy','DESCRIPTIVE_ONLY')
    require(strategy in {'DESCRIPTIVE_ONLY','EXCHANGEABILITY_REVIEWED'},'Unknown overlap-null strategy')
    if strategy=='EXCHANGEABILITY_REVIEWED':
        require(c.get('shared_gwas_dependence_justification'),'Review dependence from shared sleep GWAS/cohorts, not only shared LD blocks')
    pair_sleep=c['pair_sleep_traits'];pairs=sorted(pair_sleep)
    require(len(pairs)>=2,'Need at least two prespecified pairs')
    rows=list(read_tsv(c['units'],['pair_id','locus_id','stratum','eligible','selected','genome_build']))
    data=collections.defaultdict(dict);strata={};builds=set()
    for r in rows:
        pair=r['pair_id'];locus=r['locus_id']
        require(pair in pair_sleep and locus and r['stratum'],'Unknown pair or missing unit')
        require(locus not in data[pair],'Duplicate pair/block result')
        require(r['eligible'] in {'0','1'} and r['selected'] in {'0','1'},'Expected binary flags')
        require(r['selected']=='0' or r['eligible']=='1','Selected block must be eligible')
        require(locus not in strata or strata[locus]==r['stratum'],'Inconsistent matched stratum for same block')
        strata[locus]=r['stratum'];builds.add(r['genome_build']);data[pair][locus]=r
    require(len(builds)==1 and next(iter(builds)) in {'GRCh37','GRCh38'},'Mixed/unknown genome builds')
    output=[]
    for a,b in itertools.combinations(pairs,2):
        common={x for x in set(data[a]) & set(data[b]) if data[a][x]['eligible']==data[b][x]['eligible']=='1'}
        sa={x for x in common if data[a][x]['selected']=='1'}
        sb={x for x in common if data[b][x]['selected']=='1'}
        inter=sa & sb;union=sa | sb
        row={'pair_a':a,'pair_b':b,'same_sleep_trait':pair_sleep[a]==pair_sleep[b],
             'common_eligible_blocks':len(common),'selected_a_in_common':len(sa),'selected_b_in_common':len(sb),
             'shared_blocks':len(inter),'jaccard':len(inter)/len(union) if union else 'NA',
             'status':'TESTED','p':None,'family_fdr':None,'shared_block_ids':';'.join(sorted(inter))}
        if len(common)<int(c.get('minimum_common_blocks',20)):
            row['status']='INSUFFICIENT_COMMON_COVERAGE'
        elif not sa or not sb:row['status']='NO_SELECTED_BLOCKS_IN_ONE_OR_BOTH_PAIRS'
        elif strategy=='DESCRIPTIVE_ONLY':row['status']='DESCRIPTIVE_OVERLAP'
        else:
            groups=collections.defaultdict(set)
            for x in common:groups[strata[x]].add(x)
            ans=conditional_overlap([(len(g),len(g&sa),len(g&sb),len(g&inter)) for g in groups.values()])
            row.update(p=ans['p'],expected_overlap=ans['expected_loci'])
            if not ans['informative_strata']:row.update(status='NO_EXCHANGEABLE_CONTROLS',p=None)
        output.append(row)
    q=bh([r['p'] for r in output])
    for r,v in zip(output,q):r['family_fdr']=v
    fields=['pair_a','pair_b','same_sleep_trait','common_eligible_blocks','selected_a_in_common',
            'selected_b_in_common','shared_blocks','jaccard','expected_overlap','p','family_fdr','status','shared_block_ids']
    with transaction(root,name,stage='cross_disorder_block_comparison',inputs=[manifest,c['units']],parameters=c,
                     synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'comparisons.tsv',fields,({k:'NA' if r.get(k) is None else r[k] for k in fields} for r in output))
        meta['scientific_status']='PASS' if all(r['status'] in {'TESTED','DESCRIPTIVE_OVERLAP'} for r in output) else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'expected_contrasts':len(output),'dependence_strategy':strategy,
                  'interpretation':'Regional overlap, not evidence that causal variants or cell mechanisms are identical. Compare same-sleep and different-sleep contrasts separately.'})
    return Path(root)/name
