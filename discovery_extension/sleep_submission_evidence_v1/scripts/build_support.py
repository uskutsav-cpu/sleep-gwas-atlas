#!/usr/bin/env python3
"""Join versioned evidence, generate main tables and accurate human handoff."""
import csv
import json
from collections import Counter
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P.parents[1]

def read(path):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write(path, rows):
    if not rows: raise ValueError('empty output')
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',restval='NA'); w.writeheader();w.writerows(rows)

def main():
    t=P/'tables'; d=read(t/'global_1200.tsv'); r=read(t/'replication_family_217.tsv')
    n=read(t/'NOVELTY_MASTER_1200.tsv'); s=json.loads((t/'numerical_summary.json').read_text())
    novelty={x['pair_id']:x for x in n}
    repl={x['pair_id']:x for x in r}
    dids={x.get('pair_id',x['sleep_trait']+'__'+x['extension_trait_id']) for x in d}
    if len(novelty)!=len(n) or set(novelty)!=dids: raise ValueError('novelty universe mismatch')
    source_comp={x['pair_id']:x for x in read(t/'replication_source_comparison.tsv')}
    sleep_meta={x['trait_id']:x for x in read(t/'sleep_trait_metadata.tsv')}
    external_meta={x['extension_trait_id']:x for x in read(R/'discovery_extension/config/candidate_traits.tsv')}
    claims=[]
    for x in d:
        pid=x.get('pair_id',x['sleep_trait']+'__'+x['extension_trait_id'])
        prior=novelty[pid]; rr=repl.get(pid,{})
        sl=sleep_meta[x['sleep_trait']];ex=external_meta[x['extension_trait_id']];sc=source_comp.get(pid,{})
        claims.append({'claim_id':'GLOBAL:'+pid,'pair_id':pid,'claim_type':'OBSERVED_GLOBAL_GENETIC_CORRELATION',
             'sleep_trait':x['sleep_trait'],'extension_trait_id':x['extension_trait_id'],'rg':x['rg'],'se':x['se'],
             'p':x['p'],'bh_fdr':x['extension_fdr'],'significance':str(float(x['extension_fdr'])<.05),
             'numerical_status':'INDEPENDENT_RECALCULATION_FROM_ARCHIVED_PRECISION',
             'native_source_rerun':'BLOCKED_MISSING_GWAS_AND_LDSC_REFERENCES',
             'replication_class':rr.get('audit_replication_class','NOT_IN_SELECTED_REPLICATION_FAMILY'),
             'fully_independent_two_trait_replication':'False','novelty_assessment':prior['current_novelty_class'],
             'prior_dois':prior['prior_dois'],'current_pair_evidence_status':prior['current_pair_evidence_status'],
             'sleep_source_id':sl['source_id'],'discovery_outcome_accession':ex['study_accession'],
             'discovery_outcome_source_url':ex['source_url'],'discovery_source_checksum':ex['checksum'],
             'replication_accession':rr.get('replication_study_accession','NA'),
             'replication_source_checksum':sc.get('replication_checksum','NA'),
             'source_metadata_evidence':'tables/phenotype_and_source_metadata.tsv;tables/sleep_trait_metadata.tsv;tables/source_integrity_ledger.tsv',
             'source_comparison_evidence':'tables/replication_source_comparison.tsv#'+pid if rr else 'NOT_IN_SELECTED_REPLICATION_FAMILY',
             'discovery_source':'discovery_extension/results/ldsc/extension_rg_matrix.tsv',
             'evidence_tables':'tables/global_1200.tsv;tables/NOVELTY_MASTER_1200.tsv;tables/replication_family_217.tsv',
             'allowed_interpretation':'Global association only; no causal, local, molecular or clinical benefit inference',
             'open_objections':'Sleep reused for available validation; cohort overlap depends on pair (see source audit); priority unresolved where no direct table; original sleep h2 incomplete'})
    write(t/'full_claim_evidence_ledger.tsv',claims)
    panel=read(R/'config/analysis_panel.tsv'); ext=read(R/'discovery_extension/config/candidate_traits.tsv')
    composition=[{'analysis':'immutable_core','traits':len(panel),'sleep_traits':sum(x['domain']=='sleep' for x in panel),
                  'non_sleep_traits':sum(x['domain']!='sleep' for x in panel),'planned_tests':sum(x['domain']=='sleep' for x in panel)*sum(x['domain']!='sleep' for x in panel),
                  'verification':'HISTORICAL_OUTPUT_HASH_ONLY; MISSING_CHECKPOINT_COMPONENTS'},
                 {'analysis':'exploratory_extension','traits':len(ext)+s['sleep_traits'],'sleep_traits':s['sleep_traits'],
                  'non_sleep_traits':len(ext),'planned_tests':len(d),'verification':'NUMERICAL_RECALCULATION; NO_NATIVE_RERUN'}]
    write(t/'main_panel_composition.tsv',composition)
    metrics=[{'metric':k,'value':v,'evidence':'tables/numerical_summary.json'} for k,v in s.items() if isinstance(v,(int,float))]
    metrics.append({'metric':'candidate_phenotypes_considered','value':len(read(R/'discovery_extension/results/candidate_pool.tsv')),'evidence':'discovery_extension/results/candidate_pool.tsv'})
    metrics.append({'metric':'external_h2_primary_pass','value':len([x for x in read(R/'discovery_extension/results/ldsc/extension_trait_readiness.tsv') if x['primary_rg_eligibility']=='PRIMARY_PASS']),'evidence':'discovery_extension/results/ldsc/extension_trait_readiness.tsv'})
    for label,count in Counter(x['audit_replication_class'] for x in r).items():
        metrics.append({'metric':'replication_'+label.lower(),'value':count,'evidence':'tables/replication_family_217.tsv'})
    write(t/'main_discovery_summary.tsv',metrics)
    groups={}
    for x in d:
        z=groups.setdefault(x['phenotype_domain'],{'phenotype_domain':x['phenotype_domain'],'comparisons':0,'bh_significant':0,'distinct_phenotypes':set()})
        z['comparisons']+=1; z['bh_significant']+=float(x['extension_fdr'])<.05;z['distinct_phenotypes'].add(x['extension_trait_id'])
    for z in groups.values(): z['distinct_phenotypes']=len(z['distinct_phenotypes']);z['interpretation']='Correlated comparison counts; not independent biological discoveries'
    write(t/'main_phenotype_domains.tsv',list(groups.values()))
    protocols='''# Versioned analysis protocol

Date: 2026-10-08, client execution date. Design: retrospective evidence and numerical verification of the existing exploratory extension.

The historical 1,200-test BH family, 217-candidate Bonferroni family, source identities, outcomes and files remain immutable. Newly calculated values carry archived-precision labels. Precision perturbation bounds, estimator-covariance rho grids, displayed-intercept strata, measurement strata, phenotype matching strata and phenotype-profile similarity are exploratory audits. Profile similarity is not external phenotype genetic correlation. No new confirmatory GWAS family, conditional analysis or causal hypothesis was launched.

Native LDSC replay needs the exact 12 sleep inputs, 100 extension inputs, replication inputs, checksum-pinned HapMap3 and EUR LD/weight references and pinned LDSC runtime. Current receipts and table arithmetic cannot replace them. Full dense data/ancestry-matched signed LD would be required for selected local work. Local methods are not essential to a carefully limited global paper, but cannot rescue novelty by being added indiscriminately. No unverified LD or unrelated Brain6/frailty input is used.
'''
    (P/'ANALYSIS_PROTOCOL.md').write_text(protocols)
    (P/'14_AI_USE_AND_RESEARCH_CONTRIBUTIONS.md').write_text('''# Actual AI-use and contribution record

Codex executed repository inspection, Python/JavaScript programming, numerical verification, source and licensing audit, primary literature searches and supplement extraction, statistical criticism, and deterministic plotting. Three Codex subagents performed numerical, provenance and literature/journal work with separate file ownership. Root reviewed their outputs and incorporated cross-review objections. AI tools are not authors.

No manuscript section, Statement of Significance, cover letter or manuscript document was written or modified. No journal was contacted or submission made. No author approval, affiliation, ethics determination, funding, disclosure or human verification was invented. Figures use observed tables and standard plotting; no generative image model was used. Checklists and audit documents are technical supporting materials.

Human verification in this execution: NOT_DOCUMENTED. Model/version disclosure beyond the available Codex identity: investigator should record the app/model version visible to them. Human authors must verify the code, tables, current literature classifications, statistical interpretation and final AI-use disclosure under the current SLEEP policy.
''')
    activities=[('Repository and source inspection','Codex root and provenance reviewer'),('Statistical programming and independent implementation','Codex numerical reviewer'),
        ('Primary literature and supplementary tables','Codex literature reviewer'),('Figures and editable table package','Codex root'),('Automated adversarial criticism','Root and independent subagent cross-review'),
        ('Human scientific verification','NOT_DOCUMENTED'),('Manuscript drafting','NOT_PERFORMED')]
    write(t/'ai_use_activity_ledger.tsv',[{'activity':a,'performed_by':b,'date':'2026-10-08','human_approval':'NOT_DOCUMENTED'} for a,b in activities])
    approvals=['Author list','Author order','Corresponding author','Institutional affiliations','CRediT contributions','Funding and grant identifiers',
       'Conflicts of interest','Ethics determination for aggregate summary-statistics research','Source data-use and redistribution permissions',
       'AI-use disclosure','Repository release authorization beyond current authorized branch changes','Final scientific interpretation',
       'Independent human statistical-genetics review','Manuscript approval','Journal submission authorization']
    write(t/'human_approval_checklist.tsv',[{'required_information':a,'status':'REQUIRES_HUMAN_INPUT','authorized_person':'UNSUPPLIED','evidence':'NOT_DOCUMENTED'} for a in approvals])
    (P/'15_ETHICS_AUTHOR_APPROVAL_CHECKLIST.md').write_text('# Investigator approval requirements\n\n'+
        'All 15 required decisions are in `tables/human_approval_checklist.tsv`, with no assumed answers. Public aggregate GWAS availability does not itself establish IRB exemption, authorship, permitted redistribution or institutional approval.\n\nNo correspondence, access forms, agreements, charges, DOI minting, archive publication or submission was performed. The dedicated branch push was explicitly authorized in the user brief; further investigator-controlled releases require documented rights and author decisions.\n')
    (P/'18_MANUSCRIPT_AUTHOR_HANDOFF.md').write_text('''# Evidence map for the human author

This file maps evidence and open decisions; it contains no manuscript prose.

| Evidence | Verified package location | Author decision |
|---|---|---|
| Complete discovery universe, estimates, precision sensitivity and BH | tables/global_1200.tsv; 03_GLOBAL_RESULTS_VERIFICATION.md | Keep 1,200 denominator and exploratory scope |
| All positive and negative replication candidates | tables/replication_family_217.tsv | Retain unavailable/QC-ineligible outcomes |
| Outcome-side positives and exact/comparable definitions | tables/replicated_23.tsv; 04_REPLICATION_INDEPENDENCE_AUDIT.md | Reused sleep GWAS precludes fully independent two-trait claim |
| Heterogeneity and covariance grid | tables/HETEROGENEITY_MASTER.tsv; 05_HETEROGENEITY_VALIDATION.md | Seven nominal flags are qualified; covariance is unknown |
| Current pair-level prior art | tables/NOVELTY_MASTER_1200.tsv; PRIOR_STUDY_COMPARISON.md | Verify direct comparators; do not claim first-ever novelty |
| Source integrity and recovery | tables/source_integrity_ledger.tsv; 02_SOURCE_AND_PROVENANCE_AUDIT.md | Nine checkpoint originals remain missing; environment order drift |
| All figures and exact numerical inputs | figures/; figures/source_data/ | Select smallest effective set after human scientific review |
| Reporting and journal requirements | 11_SLEEP_AUTHOR_GUIDELINES_AUDIT.md; 12_STROBE_STREGA_CHECKLIST.md | Complete manuscript-only and investigator-only items |
| Claim boundaries | tables/full_claim_evidence_ledger.tsv | No causality, mechanism, shared variant, gene or clinical benefit claim |

Required human scientific questions: Does outcome-side validation add enough beyond Morrison 2024, Goodman 2025 and the 2026 Sleep Chart? Can exact original sleep h2 and complete frozen core outputs be recovered? Can shared-estimator covariance be estimated from verified jackknife inputs? Is a new independent sleep cohort with comparable definitions accessible and scientifically justified? Is SLEEP the appropriate target given limited measurement depth and incremental overlap?

Author/investigator actions are listed in `15_ETHICS_AUTHOR_APPROVAL_CHECKLIST.md`. Manuscript construction, significance statement, cover letter, final authorship and submission remain human tasks.
''')
    (P/'REPRODUCE.md').write_text('''# Reproduction boundary and commands

Run from the repository root using Python 3.11. The source-free numerical audit uses only the Python standard library. Figures use the exact pinned dependencies in `requirements/figures.txt`. Keep raw restricted data outside Git.

```sh
python3.11 -m venv work/evidence-env
task_python=work/evidence-env/bin/python
$task_python -m pip install -r discovery_extension/sleep_submission_evidence_v1/requirements/figures.txt -r discovery_extension/sleep_submission_evidence_v1/requirements/tooling.txt
# Verify delivered bytes before regeneration changes outputs or execution logs.
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/finalize_package.py --verify
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/verify_statistics.py
$task_python -m unittest discover -s discovery_extension/sleep_submission_evidence_v1/tests -v
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/build_figures.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/independent_statistics.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py --verify-canonical
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/build_support.py
$task_python discovery_extension/synthetic/test_extension_ldsc_collation.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/validate_package.py
```

Optional source-dependent replay is separate: `audit_provenance.py` needs the recorded local historical worktrees/Git history and acquisition receipts; `audit_literature.py --acquire` reconstructs required official-source cache bytes, followed by `audit_literature.py` with `requirements/exports.txt` dependencies. This does not repeat every current database search. Do not treat `--verify-canonical` as fresh source verification. The provenance audit supports local historical worktree paths; inspect `--help` or its script before rerunning and preserve current acquisition receipts. The literature audit uses saved extracted primary tables and retains explicit unresolved classifications; live searches must be repeated for a later execution date. The workbook builder uses `@oai/artifact-tool` from Codex bundled dependencies and reads canonical TSV tables; editable XLSX is a static research-data export, not a replacement statistical engine.

Hash verification: `shasum -a 256 -c discovery_extension/sleep_submission_evidence_v1/hashes/PACKAGE_SHA256SUMS` from the repository root. Final receipt/manifests intentionally do not hash themselves. Current literature cache acquisition is `audit_literature.py --acquire`; source-free canonical evidence checking is `audit_literature.py --verify-canonical`. These are distinct evidence levels. For the workbook, run `workbook_input.py work/workbook_data.json`, expose the bundled artifact-tool modules in the script directory, then run `build_workbook.mjs work/workbook_data.json discovery_extension/sleep_submission_evidence_v1` with bundled Node. For editable Word tables, run `build_word_tables.py` with bundled Python/python-docx, then the Documents skill's `render_docx.py` for every output before delivery.

Scientific native rerun: BLOCKED. No original sleep/extension munged GWAS or EUR LD/weight reference is local. Acquire identities and hashes from original source manifests and licenses; use the frozen versioned streaming scripts only after adequate compute/network/storage feasibility is approved. Never insert regenerated files under missing frozen artifact names. The old full dense mirror estimates about 246.91 GiB and was not initiated.

Existing historical tests with missing scientific artifacts or absent R tooling are not a native scientific pass. Initial sparse-checkout failures are retained in logs; complete-checkout baseline and new-package results are separately recorded. Historical scientific and manuscript paths are unchanged.
''')
    (P/'requirements/figures.txt').write_text('numpy==1.26.4\npandas==2.2.3\nmatplotlib==3.9.4\n')
    (P/'requirements/numerical.txt').write_text('# Python 3.11 standard library only. No third-party numerical dependency.\n')
    print(json.dumps({'claim_rows':len(claims),'main_tables':3,'approval_items':len(approvals)}))

if __name__=='__main__':main()
