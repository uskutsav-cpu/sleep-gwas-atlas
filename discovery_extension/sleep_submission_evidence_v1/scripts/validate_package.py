#!/usr/bin/env python3
"""Fail closed for package consistency, separate external scientific gates."""
import csv, hashlib, json, sys
from pathlib import Path
P=Path(__file__).resolve().parents[1]
R=P.parents[1]
REPORTS=['BASELINE_AND_REPOSITORY_AUDIT','SOURCE_AND_PROVENANCE_AUDIT','GLOBAL_RESULTS_VERIFICATION',
 'REPLICATION_INDEPENDENCE_AUDIT','HETEROGENEITY_VALIDATION','NOVELTY_AND_PRIOR_ART','SENSITIVITY_RESULTS',
 'SLEEP_PHENOTYPE_AND_CLINICAL_INTERPRETATION','NEGATIVE_AND_BLOCKED_RESULTS','ADVERSARIAL_SCIENTIFIC_REVIEW',
 'SLEEP_AUTHOR_GUIDELINES_AUDIT','STROBE_STREGA_CHECKLIST','DATA_AND_CODE_AVAILABILITY_INVENTORY',
 'AI_USE_AND_RESEARCH_CONTRIBUTIONS','ETHICS_AUTHOR_APPROVAL_CHECKLIST','PUBLIC_RELEASE_AND_LICENSE_AUDIT',
 'SLEEP_EDITORIAL_READINESS','MANUSCRIPT_AUTHOR_HANDOFF']

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def indexed(rows,key):
    m={x[key]:x for x in rows}
    if len(m)!=len(rows):raise ValueError('Duplicate '+key)
    return m

def check(P=P, require_reports=True):
    if require_reports:
        for i,name in enumerate(REPORTS,1):
            file=P/f'{i:02}_{name}.md'
            if not file.is_file() or file.stat().st_size<100:raise ValueError('Missing report '+str(file))
    d=read(P/'tables/global_1200.tsv');r=read(P/'tables/replication_family_217.tsv')
    pairs=indexed(d,'pair_id');reps=indexed(r,'pair_id')
    original=read(R/'discovery_extension/results/ldsc/extension_rg_matrix.tsv')
    expected={x['sleep_trait']+'__'+x['extension_trait_id'] for x in original}
    if set(pairs)!=expected or len(pairs)!=1200:raise ValueError('Discovery family drift')
    orig_reps=read(R/'discovery_extension/results/replication/replication_results.tsv')
    if set(reps)!={x['pair_id'] for x in orig_reps} or len(reps)!=217:raise ValueError('Replication family drift')
    for x in original:
        audited=pairs[x['sleep_trait']+'__'+x['extension_trait_id']]
        for c in x:
            if audited[c]!=x[c]:raise ValueError('Original scientific value changed '+c)
    for x in orig_reps:
        audited=reps[x['pair_id']]
        for c in x:
            if audited[c]!=x[c]:raise ValueError('Original replication value changed '+c)
        if x['replication_class']=='NO_INDEPENDENT_DATASET' and audited['replication_rg']!='NA':raise ValueError('Missing value converted')
    hits=read(P/'tables/significant_603.tsv');positive=read(P/'tables/replicated_23.tsv');heter=read(P/'tables/heterogeneity_7.tsv')
    if set(indexed(hits,'pair_id'))!={k for k,x in pairs.items() if float(x['extension_fdr'])<.05}:raise ValueError('FDR subset mismatch')
    if set(indexed(positive,'pair_id'))!={k for k,x in reps.items() if x['replication_class']=='REPLICATED'}:raise ValueError('Replication subset mismatch')
    if set(indexed(heter,'pair_id'))!={x['pair_id'] for x in positive if float(x['heterogeneity_p'])<.05}:raise ValueError('Heterogeneity subset mismatch')
    for file in ['NOVELTY_MASTER_1200.tsv','full_claim_evidence_ledger.tsv']:
        if set(indexed(read(P/'tables'/file),'pair_id'))!=expected:raise ValueError('Evidence universe drift '+file)
    if set(indexed(read(P/'tables/replication_source_comparison.tsv'),'pair_id'))!=set(reps):raise ValueError('Source comparison family drift')
    if any(x['audit_replication_class']!='EXTERNAL_OUTCOME_SIDE_REPLICATION' for x in positive):raise ValueError('Independence overclaim')
    if any(x['fully_independent_two_trait_replication']!='False' for x in read(P/'tables/full_claim_evidence_ledger.tsv')):raise ValueError('Unsupported independence')
    receipt=json.loads((P/'figures/FIGURE_EXECUTION_RECEIPT.json').read_text())
    for rel,digest in receipt['source_hashes'].items():
        if hashlib.sha256((R/rel).read_bytes()).hexdigest()!=digest:raise ValueError('Figure source hash drift')
    return {'status':'PASS_PACKAGE_CONSISTENCY','global_rows':len(d),'significant_rows':len(hits),
            'replication_family':len(r),'positive_rows':len(positive),'nominal_heterogeneity_rows':len(heter),
            'native_scientific_reproduction':'BLOCKED_EXTERNAL_INPUTS','fully_independent_two_trait_replication':0}

if __name__=='__main__':
    try:
        result=check();(P/'logs/package_consistency.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    except (ValueError,KeyError,FileNotFoundError) as error:
        print('FAIL_PACKAGE_CONSISTENCY: '+str(error),file=sys.stderr);sys.exit(1)
