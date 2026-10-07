"""Deterministic figures and investigator-review materials from saved results.

No GWAS fit is performed here. Public source tables contain derived results,
not restricted raw association subsets.
"""
from pathlib import Path
import hashlib, json, shutil, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
AREA = ROOT / 'brain6/translational_psychiatry_research_v1'
R = AREA / 'research_v1'
M = AREA / 'MANUSCRIPT'
F = M / 'figures'
F.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
    'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'brain6-research-v1',
    'axes.spines.top':False,'axes.spines.right':False})
inputs = {}
def read(path):
    inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, sep='\t')
def save(fig, name):
    fig.savefig(F / (name+'.pdf'), bbox_inches='tight', metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(F / (name+'.svg'), bbox_inches='tight', metadata={'Date':None})
    fig.savefig(F / (name+'.png'), bbox_inches='tight', dpi=180, metadata={'Software':'Brain6 saved-result figure builder'})
    plt.close(fig)

g = read(ROOT/'brain6/confirmatory_v2/MANUSCRIPT/figures/figure1_source.tsv')
g.to_csv(F/'figure1_source.tsv', sep='\t', index=False)
sleep = list(dict.fromkeys(g.sleep_trait))
disorders = ['adhd','mdd','scz','bipolar','pd','ad']
actual = list(dict.fromkeys(g.disorder))
disorders = ['adhd','mdd'] + [x for x in actual if x not in ['adhd','mdd']]
matrix = g.pivot(index='sleep_trait',columns='disorder',values='rg').reindex(index=sleep,columns=disorders)
sig = g.pivot(index='sleep_trait',columns='disorder',values='significant_fdr').reindex(index=sleep,columns=disorders)
fig, ax = plt.subplots(figsize=(8.2,6.8))
im = ax.imshow(matrix.to_numpy(), cmap='RdBu_r', vmin=-1,vmax=1, aspect='auto')
for i in range(len(sleep)):
    for j in range(len(disorders)):
        v=matrix.iloc[i,j]
        mark='*' if str(sig.iloc[i,j]).lower()=='true' else ''
        ax.text(j,i,f'{v:.2f}{mark}',ha='center',va='center',fontsize=8.5,
            color='white' if abs(v)>.55 else 'black')
labels = dict(zip(g.disorder,g.brain_disorder_label))
ax.set_xticks(range(len(disorders)),[labels[x] for x in disorders], rotation=25, ha='right')
ax.set_yticks(range(len(sleep)),[x.replace('_',' ') for x in sleep])
ax.set_title('Inherited sleep–brain disorder genetic correlations',pad=14)
fig.colorbar(im,ax=ax,shrink=.7,label='Genetic correlation')
fig.text(.02,.005,'* Original 396-test FDR significance; 72 rows shown. Detection counts are not architecture tests.',fontsize=9)
fig.tight_layout(rect=(0,.04,1,1))
save(fig,'figure1_global_map')

pips=read(R/'ld/native_finemap_variant_PIPs.tsv')
coloc=read(R/'ld/native_coloc_all_signal_pairs.tsv')
primary=pips[pips.model=='PRIMARY_L10_MEDIAN_NEFF'].copy()
primary['position_GRCh37']=primary.variant_key.str.split(':').str[1].astype(int)
primary.to_csv(F/'figure2_PIP_source.tsv',sep='\t',index=False)
coloc.to_csv(F/'figure2_coloc_source.tsv',sep='\t',index=False)
fig,axs=plt.subplots(3,1,figsize=(8.7,8.2),gridspec_kw={'height_ratios':[1,1,1.15]})
colors={'insomnia':'#2268a2','adhd':'#bd5a24'}
for ax,trait in zip(axs[:2],['insomnia','adhd']):
    q=primary[primary.trait==trait]
    ax.scatter(q.position_GRCh37/1e6,q.PIP,s=9,color=colors[trait],alpha=.8,rasterized=False)
    ax.set_xlim(103.447968,104.447968);ax.set_ylim(0,.32)
    ax.set_ylabel(f'{trait.upper()} PIP')
    ax.set_xlabel('Chromosome 5 position (GRCh37, Mb)')
axs[0].set_title('Known insomnia–ADHD region: native primary fine-mapping',pad=12)
names={'PRIMARY_L10_MEDIAN_NEFF':'Primary L=10','SENSITIVITY_L1_MEDIAN_NEFF':'L=1',
       'SENSITIVITY_L5_MEDIAN_NEFF':'L=5','SENSITIVITY_L10_N_INFINITY':'Large N'}
for model,q in coloc.groupby('model',sort=False):
    q=q.sort_values('p12')
    axs[2].plot(q.p12,q['PP.H4.abf'],'o-',label=names.get(model,model))
axs[2].set_xscale('log');axs[2].set_ylim(.8,1.005)
axs[2].set_xlabel('Shared-variant prior p₁₂');axs[2].set_ylabel('coloc-SuSiE H4')
axs[2].legend(fontsize=8,loc='lower right');axs[2].grid(alpha=.18)
fig.text(.02,.006,'Same discovery GWAS; one qualifying credible set per trait. Model support does not establish replication or mechanism.',fontsize=8.5)
fig.tight_layout(rect=(0,.04,1,1));save(fig,'figure2_native_chr5')

cal=read(R/'statistics/calibration_simulations.tsv')
cal.to_csv(F/'figure3_source.tsv',sep='\t',index=False)
fig,axs=plt.subplots(1,2,figsize=(9.2,4.7),sharey=True)
for ax,c in zip(axs,[0.,.2]):
    for h,q in cal[cal.true_intercept==c].groupby('local_h2_both'):
        q=q.sort_values('ld_ar1')
        y=q['p_alpha0.05_rate'];lo=q['p_alpha0.05_wilson_lower'];hi=q['p_alpha0.05_wilson_upper']
        ax.errorbar(q.ld_ar1,y,yerr=np.array([y-lo,hi-y]),fmt='o-',capsize=3,label=f'Local h²={h:g}')
    ax.axhline(.05,color='black',linestyle='--',linewidth=1,label='Nominal 5%')
    ax.set_title(f'Overlap intercept = {c:g}')
    ax.set_xlabel('AR(1) LD parameter');ax.set_ylim(0,.57);ax.grid(alpha=.15)
axs[0].set_ylabel('Truth-null rejection fraction');axs[1].legend(fontsize=8,loc='upper left')
fig.suptitle('Pinned-source covariance uncertainty: illustrative counterexamples',fontsize=12)
fig.text(.025,.01,'10,000 draws/scenario; exact synthetic LD, covariance zero. These are not estimated Brain6 false-positive rates.',fontsize=8.5)
fig.tight_layout(rect=(0,.06,1,.95));save(fig,'figure3_covariance_diagnostic')
(F/'figure_provenance.json').write_text(json.dumps({'classification':'SAVED_RESULT_VISUALIZATION',
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256':inputs,'numpy':np.__version__,'pandas':pd.__version__,
    'matplotlib':matplotlib.__version__,'deterministic_metadata':True},indent=2)+'\n')

# Preserve complete derived source tables needed for the inherited manuscript.
ST=M/'source_tables';ST.mkdir(exist_ok=True)
for path in (ROOT/'brain6/confirmatory_v2/MANUSCRIPT/source_tables').glob('*'):
    if path.is_file():shutil.copy2(path,ST/path.name)
for src,name in [(AREA/'06_REPLICATION_RESULTS.tsv','clinical_fixed_variant_slots.tsv'),
    (R/'ld/native_finemap_fit_summary.tsv','native_finemap_fit_summary.tsv'),
    (R/'ld/native_finemap_credible_sets.tsv','native_credible_sets.tsv'),
    (R/'ld/native_coloc_all_signal_pairs.tsv','native_coloc_all_signal_pairs.tsv'),
    (R/'global_contrasts/contrasts.tsv','global_contrasts.tsv'),
    (R/'MOLECULAR_EVIDENCE_MASTER.tsv','molecular_evidence.tsv')]:shutil.copy2(src,ST/name)

abstract='''Cross-trait association does not establish shared causal variants or molecular mechanisms. We reassessed a sleep–brain disorder genetic atlas using source provenance, public genotype-derived linkage disequilibrium (LD), native multiple-signal fine-mapping, clinical phenotype sensitivity and explicit statistical diagnostics. The inherited map contained 72 comparisons and 35 original-family significant entries; 25 screening candidates occupied 20 regions. Canonical LAVA remained failed, with 3,720 unestimated of 17,465 cells against its frozen ceiling of 873. Three newly reconstructed European reference-LD matrices passed separate validation gates. At a previously reported chr5 insomnia–ADHD region, 2,185 shared variants supported convergent SuSiE fits and one qualifying credible set per trait. Coloc-SuSiE H4 was 0.9916 at the default prior and 0.9216 at the lower prior, with unresolved causal-variant identity. Two of four fixed-variant clinical pair slots passed the retained family; one did not and one was unestimated. These checks remained phenotype sensitivities with unverified participant independence. Accessible prior supplements supported regional evidence for 22 protected candidates. Truth-known diagnostics identified a source-specific mismatch between covariance estimation and reported uncertainty, without estimating empirical false-positive rates. No molecular effector was established. The results strengthen a known shared-association model while retaining historical QC failures and unresolved inference limits; they do not establish a new mechanistic discovery or independent phenotype-matched replication.'''
assert 150<=len(abstract.split())<=250
(M/'MP_ABSTRACT.txt').write_text(abstract+'\n')
declarations=(ROOT/'brain6/confirmatory_v2/MANUSCRIPT/DECLARATIONS_AND_AUTHOR_FIELDS.md').read_text()
declarations+='\nThe new experiment recovered exact insomnia/ADHD sources and executed native analyses. Full historical reproduction still requires original LAVA/overlap provenance; successful new native fits do not close that hold. ADHD source-specific no-reposting terms take precedence over generic Figshare metadata. New filtered association slices, native input vectors and RDS objects remain local-only. Human permission verification is required before external publication.\n\nAI-use update: three parallel automated agents assisted reference reconstruction, source/novelty research and independent numerical/statistical review. Their checks are computational cross-checks and are not human peer review.\n'
(M/'DECLARATIONS_AND_AUTHOR_FIELDS.md').write_text(declarations)
(M/'JOURNAL_SUBMISSION_CHECKLIST.md').write_text('''# Journal-specific submission preparation

Checked 7 October 2026. Genome Medicine requires a structured abstract of no more than 350 words, three to ten keywords, research sections and complete declarations. This source supplies those sections and six keywords. Its [Mental health and neuropsychiatric disorders collection](https://link.springer.com/collections/dbjjedaefg) has a 12 May 2027 deadline; topical relevance does not meet its contribution threshold. [Research requirements](https://link.springer.com/journal/13073/submission-guidelines/research).

[Molecular Psychiatry's current HTML guidance](https://www.nature.com/mp/authors-and-referees/preparation-of-articles) lists 5,000 article words, six display items, 100 references and a 150–250-word unstructured abstract. Its older linked PDF differs. This package uses conservative 3,500 main words, four display items (three figures and one table) and 15 references. MP_ABSTRACT.txt is the unstructured alternative. Investigators must confirm current limits on submission day.

Genome Medicine / Molecular Psychiatry original-discovery verdict: NO_GO. Specialist fallback BMC Medical Genomics: scope is plausible for a bounded secondary genomic study, but human statistical review and editorial assessment are required. No publication probability is asserted.

Human submission approval is held on: original-source audit closure; covariance inference review; exact data permissions; responsible authors, funding, ethics and conflicts; all author approvals; permanent archive/version selection; journal formatting and final scientific judgment. Do not remove the investigator-review notice until those are complete. No submission or third-party communication occurred.
''')
(M/'COVER_LETTER_DRAFT.md').write_text('''# Draft cover letter — investigator completion; not sent

Dear Editor,

Please consider “A source-qualified reassessment of locus sharing in a sleep–brain disorder genetic atlas” as a secondary genomic validation study. The work preserves an inherited 72-comparison atlas, 25 candidate associations and failed canonical local-analysis gate. A new public-genotype reference permits native multiple-signal analysis of a known insomnia–ADHD region, while fixed-variant clinical sensitivity and source-specific statistical diagnostics delimit its interpretation. All negative, unestimated and ineligible outcomes are retained.

The principal contribution is an auditable account of evidence progression and inference limits. We do not claim a novel causal locus, effector gene, independent phenotype-matched replication, clinical utility or a resolved correction to covariance inference. The automated checks are not independent human peer review. Genome Medicine and Molecular Psychiatry are presently NO_GO for the original-discovery claim; the investigators should select an appropriate specialist venue after scientific review.

Responsible authors, affiliations, corresponding author, permissions, ethics, funding, competing interests, originality/exclusive-submission declarations, previous dissemination and all author approvals require investigator completion. No certification is supplied by this draft. Permanent data/code archive: investigator selection required. AI assistance is disclosed in the manuscript.

Sincerely,
Corresponding author: human completion required.
''')
(M/'SUPPLEMENTARY_METHODS_AND_RESULTS.md').write_text('''# Supplementary Methods, Results and provenance

S1. Historical baseline and sources. The original 72-row display retains the original 396-test FDR family; 25 candidates span 20 display regions. Source tables accompany this directory. Historical tracked hashes bind 5,666 files. No historical threshold or result was overwritten. `../research_v1/BASELINE_AND_RESOURCES.md` and `../SOURCE_AND_ACCESS_LEDGER.tsv` supply access and inventory details.

S2. Canonical local-validation failure. LAVA has 13,745 TESTED and 3,720 NOT_RUN among 17,465 cells, above its frozen 873 ceiling. The original 88-locus pilot was not promoted. Failed source/heritability eligibility is not biological absence. Secondary SUPERGNOVA has 3,304 numerical estimates within the 8,465-slot family, three recorded positives, and zero protected same-pair overlaps. Calibration is uncertified.

S3. New empirical reference and native analysis. `../research_v1/LD_RECONSTRUCTION_REPORT.md`, `LD_VALIDATION.tsv`, and `ld/REPRODUCE.md` describe three genotype-derived matrices, 503 European samples, variant/reference checks and exact commands. `ld/FINEMAP_PRE_FIT_PROTOCOL.md` binds the new chr5 experiment. Every PIP, credible set, convergence/ELBO record, 12 signal-pair/prior results, three separate ABF sensitivities and all warnings are retained. Native inputs and RDS objects are local-only under source terms, with hashes and acquisition instructions. Independent RDS/configuration review is in `statistics/INDEPENDENT_NATIVE_FINE_MAPPING_REVIEW.md`.

S4. Clinical sensitivity. Four fixed variant/pair slots use maximum native P and Bonferroni four. Two pass, one does not and one is unestimated. Native endpoint fields, nulls, correction/amendment, exact keys, source hashes and independent raw-JSON arithmetic checks are in `sources/`. These are phenotype sensitivities, with independence unverified and no complete regional colocalization data. Exact unsent investigator requests are included.

S5. Covariance diagnostics. `statistics/calibration_protocol.json` precedes 24 × 10,000 truth-known replicates. Unchanged upstream fixtures, full historical estimator/source audit, Gaussian fourth moments, heavy-tail diagnostics and 200,000 separate variance-review draws are retained. No synthetic factor corrects historical SE/P. All 3,304 historical estimates, out-of-bounds derived correlations and missingness remain visible. See `statistics/SUPERGNOVA_REASSESSMENT.md` and `STATISTICAL_ADVERSARIAL_REVIEW.md`.

S6. Contrast and molecular sensitivity. The 72 contrasts and 648 covariance-grid rows are source-defined, retrospective and conditional on marginal calibration. Sixteen conservative contrasts pass Bonferroni 72. Four molecular components have 60 prior/scale conditions and a separate 60-condition nonpalindromic stress grid; no H4 reaches 0.8. Sixty source-function R concordance checks do not constitute a fresh molecular dataset. Candidate-C historical QTL contexts fail admission. Old 20-region enrichment failed controls; separate 18-region sensitivity had no significant tissue/pathway. No gene, cell or pathway mechanism is inferred.

S7. Novelty and limitations. Five fully paginated focused searches yielded 488 unique records and 23 prioritized competitors, with accessible supplements inspected. Regional evidence exists for 22 of 25 candidates; three Parkinson regions remain unresolved, not novel by default. Updated chr11 regional/weak-colocalization evidence and chr5 high-LD tagging are recorded in `NOVELTY_MASTER.tsv` and `COMPETING_STUDIES.md`. Source allele/access and participant-linkage gaps remain.

S8. Reproducibility and reporting. `../REPRODUCE.md`, package checksums, source/code receipts, software locks, source-free tests and separately reported local empirical integration checks distinguish engineering from biological validation. STROBE/STREGA items are supplied. Human statistical, author, ethics and source-permission review is outstanding; this is an investigator-review package.
''')
check=pd.read_csv(ROOT/'brain6/confirmatory_v2/MANUSCRIPT/STROBE_STREGA_CHECKLIST.tsv',sep='\t')
check['manuscript_location']=check.manuscript_location.str.replace('Figures1–4','Figures 1–3',regex=False).str.replace('Figure2','Supplement S2',regex=False)
check.loc[check.item==16,'status_or_limitation']='Inherited source effects; native PIPs/CS/H3–H4; clinical fields and all held outcomes; no risk prediction'
check.loc[check.item==17,'status_or_limitation']='Locked prior/model, molecular, truth-known calibration and conditional contrast sensitivities remain separate'
check.to_csv(M/'STROBE_STREGA_CHECKLIST.tsv',sep='\t',index=False)
shutil.copy2(ROOT/'brain6/confirmatory_v2/MANUSCRIPT/STREGA_ADDITIONAL_ITEMS.md',M/'STREGA_ADDITIONAL_ITEMS.md')
refs=pd.read_csv(ROOT/'brain6/confirmatory_v2/MANUSCRIPT/REFERENCES_VERIFIED.tsv',sep='\t')
refs.loc[refs.key=='schipper','title']='Shared effector genes of insomnia, anxiety and depression implicate synaptic processes as transdiagnostic drug targets'
refs.loc[refs.key=='schipper','verification']='OFFICIAL_MEDRXIV_V2_JATS; 2026-03-11; PREPRINT_NOT_PEER_REVIEWED'
for key,doi,title,src in [('rss','10.1371/journal.pgen.1010299','Fine-mapping from summary data with the Sum of Single Effects model','PLOS_GENETICS_PRIMARY_PUBLISHER'),('finngen','10.1038/s41586-022-05473-8','FinnGen provides genetic insights from a well-phenotyped isolated population','NATURE_PRIMARY_PUBLISHER')]:
    refs.loc[len(refs)]={'key':key,'doi':doi,'title':title,'metadata_id':'PRIMARY_PUBLISHER','metadata_source':src,'verification':'PRIMARY_PAGE_CHECKED_2026-10-07'}
refs.to_csv(M/'REFERENCES_VERIFIED.tsv',sep='\t',index=False)
print(json.dumps({'figures':3,'abstract_words':len(abstract.split()),'references':len(refs),'figure_inputs':len(inputs)}))
