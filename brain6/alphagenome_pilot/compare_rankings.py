#!/usr/bin/env python3
"""Descriptive within-region AVI vs single-trait ABF rank comparison.

Not a trained model, statistical test, shared-PIP, or validation analysis.
"""
import collections,csv,hashlib,json,math
from pathlib import Path
P=Path(__file__).resolve().parent/'results'; A=P/'avi_mirror'; SOURCE=A/'variant_evidence_avi.tsv'; REGIONS=P/'regions.tsv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def ranks(items, key):
    # rank-min with exact ties; list order never breaks ties.
    vals={item['snp']:key(item) for item in items if key(item) is not None}
    return {s:1+sum(other>v for other in vals.values()) for s,v in vals.items()}
def parse(x):
    if not x:return None
    y=float(x)
    if not math.isfinite(y):raise ValueError('Non-finite value')
    return y
rows=read(SOURCE); regions=read(REGIONS); groups=collections.defaultdict(list)
for row in rows:
    if row['avi_status']!='SCORED_MIRROR_EXPLORATORY':raise ValueError('Unexpected AVI status')
    groups[row['candidate_locus_id']].append(row)
rankrows=[]; summaries=[]
for region in regions:
    cid=region['candidate_locus_id']; group=groups[cid]; ar=ranks(group,lambda x:int(x['avi_raw_e5'])); trait_ranks={}
    for trait in ('insomnia','adhd'):
        valid=region[trait+'_set_status']=='EXPLORATORY_ABF_CREDIBLE_SET'
        # All variants in the union are considered, but only if this trait's credible-set gate passed.
        trait_ranks[trait]=ranks(group,lambda x:parse(x[trait+'_conditional_abf_pip'])) if valid else {}
    for row in group:
        rankrows.append({'candidate_locus_id':cid,'snp':row['snp'],'atlas_variant':row['atlas_variant'],'atlas_url':row['atlas_url'],'avi_raw_e5':row['avi_raw_e5'],'avi_phred_lower':row['avi_phred_lower'],'avi_phred_upper':row['avi_phred_upper'],'avi_rank_within_union':ar[row['snp']], 'insomnia_rank_within_union':trait_ranks['insomnia'].get(row['snp'],''),'adhd_rank_within_union':trait_ranks['adhd'].get(row['snp'],''),'insomnia_conditional_abf_pip':row['insomnia_conditional_abf_pip'],'adhd_conditional_abf_pip':row['adhd_conditional_abf_pip'],'insomnia_credible_set_member':row['insomnia_credible_set_member'],'adhd_credible_set_member':row['adhd_credible_set_member']})
    for trait in ('insomnia','adhd'):
        tr=trait_ranks[trait]; limit=min(5,len(group))
        # Inclusive rank boundary retains all ties; denominators disclose any tie expansion.
        atop={s for s,r in ar.items() if r<=limit}; ttop={s for s,r in tr.items() if r<=limit}
        summaries.append({'candidate_locus_id':cid,'trait':trait,'trait_gate':region[trait+'_set_status'],'n_union':len(group),'n_with_trait_pip':len(tr),'rank_boundary':limit,'n_avi_selected':len(atop),'n_trait_selected':len(ttop) if tr else '', 'n_overlap':len(atop&ttop) if tr else '', 'overlap_snps':';'.join(sorted(atop&ttop)) if tr else '', 'analysis_label':'DESCRIPTIVE_RANK_OVERLAP_ON_SELECTED_CREDIBLE_SET_UNION' if tr else 'NOT_ESTIMATED_TRAIT_GATE'})
for filename,data in [('variant_ranks.tsv',rankrows),('ranking_overlap.tsv',summaries)]:
    with (A/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
case=[x for x in summaries if 'chr5_103447968_104447968' in x['candidate_locus_id']]
lead=next(x for x in rankrows if x['snp']=='rs2431108' and 'chr5_103447968_104447968' in x['candidate_locus_id'])
lines=['# AVI and GWAS ranking comparison','','This descriptive comparison asks whether overall AVI and the existing trait-specific ABF PIPs prioritize the same selected variants. The comparison is restricted to each region’s union of reported 95% credible sets; it is not genome-wide and cannot validate shared causality. No machine-learning model was trained.','','Within each union, variants are ranked separately by raw AVI and conditional ABF PIP. The comparison uses a rank boundary of five (or the full set if smaller); ties share the same minimum rank and can increase the number selected. Trait rankings are **not estimated** when the original per-trait association-posterior gate did not support a credible set. Empty PIPs remain missing.','',f"At the chromosome 5 case-study locus, rs2431108 ranks {lead['avi_rank_within_union']} by AVI, {lead['insomnia_rank_within_union']} by insomnia ABF PIP, and {lead['adhd_rank_within_union']} by ADHD ABF PIP among the 71 selected variants. A shared high ranking is a hypothesis-prioritization cue, not proof of a shared causal variant.",'','| Region, GRCh37 | Trait | Candidate variants | AVI at boundary | Trait at boundary | Overlap | Variants shared |','|---|---|---:|---:|---:|---:|---|']
for row in summaries:
    cid=row['candidate_locus_id'].replace('placo_insomnia__adhd_','').replace('_',':',1).replace('_','–')
    if row['analysis_label']=='NOT_ESTIMATED_TRAIT_GATE': vals=['not estimated','not estimated','not estimated']
    else:vals=[str(row['n_trait_selected']),str(row['n_overlap']),row['overlap_snps'] or 'none']
    lines.append(f"| {cid} | {row['trait']} | {row['n_union']} | {row['n_avi_selected']} | {vals[0]} | {vals[1]} | {vals[2]} |")
lines += ['', 'These selected sets vary in size and composition. PIP ranks are conditioned on a single-causal-variant ABF model, and the AVI mirror is independently spot-checked but not verified record by record. Apparent overlap can reflect the initial selection, linkage disequilibrium or chance. The numbers have no multiple-testing-adjusted significance, and the original Brain6 evidence classes remain unchanged.','']
(A/'RANK_COMPARISON.md').write_text('\n'.join(lines))
receipt={'evidence_sha256':sha(SOURCE),'region_sha256':sha(REGIONS),'outputs_sha256':{n:sha(A/n) for n in ('variant_ranks.tsv','ranking_overlap.tsv','RANK_COMPARISON.md')},'ranking':'within candidate credible-set union; descending raw AVI or conditional ABF PIP; min rank for ties; boundary min(5, N); only original credible-set-positive traits','interpretation':'descriptive; not training, validation, inference, or shared-causal probability'}
(A/'ranking_provenance.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print('Compared',len(groups),'regions;',len(rankrows),'variant rows;',len(summaries),'region-trait states')
