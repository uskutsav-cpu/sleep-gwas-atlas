"""Build current source, component and claim ledgers without changing historical files."""
from pathlib import Path
import hashlib,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
A=ROOT/'brain6/translational_psychiatry_research_v1';R=A/'research_v1'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old=ROOT/'brain6/confirmatory_v2/SOURCE_AND_ACCESS_LEDGER.tsv'
if not old.exists():old=A/'04_DATA_AND_COHORT_LEDGER.tsv'
d=pd.read_csv(old,sep='\t',dtype=str).fillna('')
d['provenance_scope']='HISTORICAL_OR_PRIOR_FEASIBILITY; current updates explicit below'
receipts=json.loads((R/'discovery/acquisition_receipts.json').read_text())
for receipt in receipts:
    t=receipt['trait'];mask=d.trait.eq(t)
    d.loc[mask,'download_date']='2026-10-07'
    d.loc[mask,'current_readiness']='EXACT_RAW_ARCHIVE_AND_REGIONAL_SLICES_VERIFIED; original full harmonized pipeline not replayed'
    d.loc[mask,'continuation_access']='CURRENT_LOCAL_SOURCE_ACQUIRED_UNDER_SOURCE_TERMS'
    d.loc[mask,'current_hash_verification']='CURRENT_RAW_REHASH_MATCH; '+receipt['source_sha256']
    d.loc[mask,'raw_file_path']='work/research-discovery-raw/'+t+'.txt.gz'
    d.loc[mask,'source_archive_path']='work/research-discovery-raw/'+t+'.txt.gz'
    d.loc[mask,'raw_file_SHA256']=receipt['source_sha256']
    d.loc[mask,'source_archive_status']='VERIFIED_EXACT_ARCHIVE_2026-10-07'
    d.loc[mask,'provenance_scope']='CURRENT_SOURCE_RECOVERY; historical post-QC count retained separately'
    d.loc[mask,'raw_source_variant_count']=str(receipt['rows_scanned'])
    d.loc[mask,'INFO']='PRESENT_IN_RECOVERED_SOURCE; required INFO >=0.9 for new chr5 analysis'
    d.loc[mask,'per_variant_N']='SOURCE_N_PRESENT; exact new N semantics in ld/FINEMAP_PRE_FIT_PROTOCOL.md'
    d.loc[mask,'effect_metric']='log(OR); source SE is log-OR standard error'
    d.loc[mask,'effect_allele_convention']='A1; new chr5 inputs oriented to verified reference ALT'
    d.loc[mask,'EAF']='Insomnia source MAF interpreted via source minor-allele convention; ADHD control FRQ_U_186843; see native harmonization diagnostics'
    d.loc[mask,'FinnGen_overlap']='UNVERIFIED; see sources/SOURCE_COHORT_GRAPH.tsv'
    if t=='adhd':d.loc[mask,'license']='Source README forbids public reposting/secondary distribution of GWAS result files; raw regional/native inputs local-only; human publication permissions review outstanding'

feas=pd.read_csv(R/'sources/SOURCE_FEASIBILITY.tsv',sep='\t',dtype=str).fillna('')
new=[]
for row in feas.to_dict('records'):
    new.append({'trait':row['trait']+'_current_source_feasibility','study':row['source_id'],
        'DOI':row['doi'],'GWAS_accession':row['source_accession'],'source_URL':row['source_url'],
        'source_page_URL':row['phenotype_url'],'download_date':'2026-10-07',
        'sample_size':row['n_total'],'cases':row['cases'],'controls':row['controls'],
        'ancestry':row['ancestry'],'age':row['age'],'sex_composition':row['sex'],
        'cohorts':row['cohorts'],'genome_build':row['build'],'effect_allele_convention':row['effect_allele'],
        'effect_metric':row['effect_metric'],'variant_count':row['coverage'],'per_variant_N':row['variant_n'],
        'license':row['terms'],'current_readiness':row['eligibility'],'replication_class':row['classification'],
        'FinnGen_overlap':row['independence'],'current_hash_verification':row['source_sha256'],
        'release':row['release'],'phenotype_definition':row['phenotype_definition'],
        'phenotype_sha256':row['phenotype_sha256'],'effective_N_count_formula':row['neff_count_formula'],
        'N_semantics':row['n_semantics'],'terms_URL':row['terms_url'],
        'provenance_scope':'CURRENT_PUBLIC_SOURCE_AND_ACCESS_REVIEW'})
for c in [5,6,11]:
    new.append({'trait':'REFERENCE_ONLY','study':f'1000G_phase3_v5b_chr{c}_EUR',
        'release':'20130502 v5b; current official public byte ranges 2026-10-07',
        'source_URL':f'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',
        'download_date':'2026-10-07','sample_size':'503','ancestry':'EUR CEU/FIN/GBR/IBS/TSI',
        'genome_build':'GRCh37 hs37d5/b37','effect_allele_convention':'Exact REF/ALT; ALT 0/1/2 dosage',
        'variant_count':str(len(pd.read_csv(R/f'ld/chr{c}_variants.tsv',sep='\t'))),
        'license':'Official public 1000 Genomes reference; acknowledgement/terms in LD_SOURCE_MANIFEST.json',
        'current_readiness':'COMPLETE_VALIDATED_REFERENCE; no trait-GWAS independence inference',
        'replication_class':'INELIGIBLE','current_hash_verification':sha(R/f'ld/chr{c}_genotype_factors.npz'),
        'provenance_scope':'CURRENT_EMPIRICAL_REFERENCE; all source-range hashes in LD_SOURCE_MANIFEST.json'})
qtl=ROOT/'brain6/results/brain6_exploratory_functional_v1/source_gate.tsv'
for row in pd.read_csv(qtl,sep='\t').query("source_kind == 'GTEX_QTL_REGIONAL'").to_dict('records'):
    new.append({'trait':row['phenotype_or_tissue'],'study':'GTEx_v8_eQTL_Catalogue_'+row['source_id'],
        'GWAS_accession':row['source_id'],'release':'Historical GTEx v8; eQTL Catalogue QTS000015; release7 metadata',
        'source_URL':row['url'],'genome_build':'Native GRCh38; historical exact GRCh37-to-GRCh38/allele maps',
        'effect_metric':'QTL Catalogue linear-regression beta/SE; tissue-scale assumptions retained',
        'per_variant_N':row['n_semantics'],'variant_count':str(row['local_rows'])+' historic regional rows',
        'license':'Public regional tabix access was used historically; original consent/source reuse permissions require investigator review; do not republish combined ADHD association inputs',
        'current_readiness':'HISTORICAL_INPUT_SENSITIVITY; no new independent QTL acquisition',
        'current_hash_verification':sha(qtl)+' source-gate table; individual slices bound in molecular/sensitivity_protocol.json',
        'replication_class':'INELIGIBLE','provenance_scope':'HISTORICAL_MOLECULAR_SOURCE_WITH_CURRENT_NUMERICAL_SENSITIVITY'})
d=pd.concat([d,pd.DataFrame(new)],ignore_index=True).fillna('')
d.to_csv(A/'SOURCE_AND_ACCESS_LEDGER.tsv',sep='\t',index=False)
d.to_csv(A/'04_DATA_AND_COHORT_LEDGER.tsv',sep='\t',index=False)

gate=pd.read_csv(ROOT/'brain6/results/brain6_exploratory_functional_v1/coloc_gate.tsv',sep='\t')
c=gate[gate.region_grch37.str.startswith('chr5:103447968-')].copy()
assert len(c)==28 and not (pd.to_numeric(c.min_qtl_p_shared,errors='coerce')<=5e-8).any()
c['new_use']='HISTORICAL_ELIGIBILITY_ONLY; broader old window encloses new native C window; no new source experiment'
c.to_csv(R/'molecular/candidate_C_historical_QTL_admission.tsv',sep='\t',index=False)
(R/'FUNCTIONAL_EVIDENCE_REPORT.md').write_text('''# Functional evidence and advancement decision

Four retained ADHD-component molecular inputs replay their original H0–H4 vectors within 1e-10. A locked five-prior by three-QTL-scale grid executed 60 conditions; maximum H4 was 0.6175746925 and none reached 0.8. A separately frozen nonpalindromic exclusion stress grid retained all four contexts and 60 conditions; maximum H4 was 0.45525975 and none reached 0.8. Sixty independent source-function R checks verify posterior arithmetic. These are historic-input sensitivities, not independent molecular replication. Three gene IDs are represented; none is an admitted causal effector.

The new native chr5 shared-association model permits an eligibility assessment, but not automatic gene promotion. The exact 28 historical QTL context rows covering the broader chr5 interval are copied into `molecular/candidate_C_historical_QTL_admission.tsv`. None has shared QTL P <=5e-8 (best 7.82503e-6); low-count and empty splicing contexts remain visible. The broader old window encloses the new native window, so restricting it cannot create a stronger source P. These resources fail strong-QTL/coverage admission and are not rerun as a molecular headline.

The retained original inputs contain 251/1624 or 236/1505 palindromic variants. ABF sign invariance does not verify variant/strand identity. Original inputs/outcomes remain immutable; newly filtered association slices stay local-only under ADHD terms.

No effector, cell type, pathway, therapeutic target or clinical predictor is established. Historical 20-region enrichment failed its control gate. The separate 18-region positional sensitivity had no FDR-positive tissue/pathway. No new enrichment is justified without an admitted target set, matched controls and a tested background. Expression, proximity and regulatory overlap remain annotations. Dense eligible brain-cell/regulatory QTL and independent functional validation require outside data or investigator work; exact limitations are in the molecular/gene and blocker ledgers.
''')

components=[
('Historical global map','COMPLETE_VALIDATED','72 rows;35 original-family significant','MANUSCRIPT/source_tables/BRAIN6_FINAL_GLOBAL.tsv','Inherited calibrated-model estimates; no causal inference'),
('Protected PLACO candidates','COMPLETE_EXPLORATORY','25 pair candidates/20 display regions','MANUSCRIPT/source_tables/BRAIN6_FINAL_CANDIDATES.tsv','Screened same-discovery candidates; not confirmed causal loci'),
('Canonical LAVA','FAILED_QC','13745 TESTED;3720 NOT_RUN>873;17465 total','MANUSCRIPT/source_tables/BRAIN6_FINAL_LAVA.tsv','Immutable FAILED_QC_NOT_PROMOTED'),
('New full confirmatory LAVA','BLOCKED_DATA','No new admitted native family','13_RESOURCE_AND_ACCESS_BLOCKERS.md','Original source/N/overlap semantics unresolved'),
('Historical source-N pilot','COMPLETE_EXPLORATORY','88 loci;advancement failed','01_SCIENTIFIC_BASELINE.md','Not spliced into original family'),
('Historical SUPERGNOVA estimates','NEEDS_REVIEW','3304 estimates;3 recorded family positives;0 protected overlaps','research_v1/statistics/SUPERGNOVA_REASSESSMENT.md','Calibration uncertified; no empirical correction applied'),
('Three long-sleep secondary pairs','NOT_APPLICABLE','3 prospective inapplicable slots','01_SCIENTIFIC_BASELINE.md','Source/model ambiguity'),
('Original SuSiE LD family','FAILED_QC','0/32 original LD gates pass','07_LD_AND_LOCAL_VALIDATION.md','Historical decision preserved'),
('New public genotype LD','COMPLETE_VALIDATED','3 references;503 EUR;300 pair checks;3 exact rebuilds','research_v1/LD_VALIDATION.tsv','Reference validation only; not independent trait replication'),
('Exact insomnia/ADHD source recovery','COMPLETE_VALIDATED','2 archives;6 slices;exact frozen source hashes','research_v1/discovery/acquisition_receipts.json','Raw subsets local-only; original harmonized full pipeline not replayed'),
('Native candidate C SuSiE','COMPLETE_EXPLORATORY','8 converged fits;2185 SNPs;1 qualifying CS/trait','research_v1/ld/native_finemap_fit_summary.tsv','Known region; approximate single-N binary RSS'),
('Native candidate C coloc-SuSiE','COMPLETE_EXPLORATORY','12 signal-pair/prior conditions;3 ABF sensitivities','research_v1/ld/native_coloc_all_signal_pairs.tsv','H4 default .991568/lower .921630; same GWAS'),
('Candidate A/B analogous fine-mapping','NOT_RUN','Insomnia min P fails predeclared strong-association gate','research_v1/ld/FINEMAP_RESULTS_AND_LIMITATIONS.md','Gate failure is not biological null'),
('FinnGen clinical fixed-point family','COMPLETE_EXPLORATORY','4 slots;2 pass;1 nonpass;1 not estimated','06_REPLICATION_RESULTS.tsv','PHENOTYPE_SENSITIVITY; independence UNVERIFIED'),
('Independent phenotype-matched two-trait replication','BLOCKED_ACCESS','0 certified','research_v1/REPLICATION_MASTER.tsv','Need exact sources and participant/control linkage'),
('Covariance method counterexamples','COMPLETE_VALIDATED','24x10000 truth-known draws;12 unchanged-source fixtures','research_v1/statistics/calibration_simulations.tsv','Synthetic inference counterexample, not empirical false-positive estimate'),
('Independent analytic variance check','COMPLETE_VALIDATED','Dense/eigen identity;2x100000 independent MC draws','research_v1/independent_statistics/independent_variance_results.json','Checks exact synthetic formula only'),
('Covariance correction/full empirical impact','NOT_RUN','No certified correction or refitted full empirical family','12_ADVERSARIAL_REVIEW.md','Requires statistical collaborator and justified new benchmark'),
('Global source-defined contrasts','COMPLETE_EXPLORATORY','72 contrasts;16 conditional bounded differences;648 grid rows','research_v1/global_contrasts/contrasts.tsv','Same-source retrospective; not causal architectural mechanism'),
('Molecular posterior sensitivity','COMPLETE_EXPLORATORY','120 conditions;60 R arithmetic checks;0 H4>=.8','research_v1/MOLECULAR_EVIDENCE_MASTER.tsv','Historic single-component/scale assumptions'),
('Candidate C molecular admission','NOT_RUN','28 contexts;0 strong-QTL admitted','research_v1/molecular/candidate_C_historical_QTL_admission.tsv','Existing source admission fails; no new molecular null'),
('Historical 20-region enrichment','FAILED_QC','Controls4/7 vsrequired10','09_FUNCTIONAL_GENOMICS.tsv','Not promoted'),
('Historical 18-region enrichment sensitivity','COMPLETE_EXPLORATORY','0 tissue/pathway FDR positives','research_v1/CELLTYPE_AND_PATHWAY_VALIDATION.tsv','Separate positional subset; no mechanism'),
('Novelty/source primary audit','COMPLETE_VALIDATED','488 unique search records;22/25 regional priors;27 build checks;14 source-accounting checks','research_v1/NOVELTY_MASTER.tsv','Focused search not exhaustive proof;3 PD unresolved'),
('Xue main/supplement access','BLOCKED_ACCESS','Full source unavailable;metadata reviewed','research_v1/UNRESOLVED_LITERATURE_ACCESS.tsv','Do not infer novelty from access gap'),
('Human statistical and author review','NEEDS_REVIEW','No human approval asserted','MANUSCRIPT/DECLARATIONS_AND_AUTHOR_FIELDS.md','Automated sister-agent review is not human peer review'),
('Manuscript/figures/submission materials','COMPLETE_VALIDATED','Editable TeX;compiled PDF;3 vector figures;15 references;205-word MP abstract','MANUSCRIPT/brain6_reassessment.tex','Technical completion only;human scientific/submission hold'),
('Restricted/protected release, contacts, merge, submission','NOT_APPLICABLE','No actions performed','FINAL_HANDOFF.md','Outside authorized scope')]
pd.DataFrame(components,columns=['component','status','quantitative_evidence','artifact','limit']).to_csv(A/'COMPONENT_INVENTORY.tsv',sep='\t',index=False)

claims=[
('C01','72 comparisons / 35 inherited significant','MANUSCRIPT/source_tables/BRAIN6_FINAL_GLOBAL.tsv','Original FDR-family display','Inherited original-family significance; no new FDR','Historical marginal/source assumptions','No independent replication'),
('C02','25 pair candidates / 20 regions','MANUSCRIPT/source_tables/BRAIN6_FINAL_CANDIDATES.tsv','Protected PLACO screen','Candidate associations, not confirmed causal loci','Post-screen region rule;source overlap','Same discovery GWAS'),
('C03','3720 NOT_RUN exceeds 873','MANUSCRIPT/source_tables/BRAIN6_FINAL_LAVA.tsv','Frozen family eligibility','Canonical LAVA FAILED_QC_NOT_PROMOTED','Source/model eligibility failure','Immutable historic gate'),
('C04','3 valid new LD references','research_v1/LD_VALIDATION.tsv','Exact genotype Pearson LD and numeric gates','New references satisfy their separate identity/numeric protocol','503 EUR samples;rank502;build and allele contracts','300 separate pairs and exact factor rebuilds'),
('C05','8 fits;1 CS per trait;2185 shared variants','research_v1/ld/native_finemap_fit_summary.tsv','Native susieR0.14.2','Reference-qualified exploratory known-region support','Single-N binary approximation;reference/sources;diffuse PIPs','Independent RDS/PIP/CS/ELBO recomputation, not new cohorts'),
('C06','H4 default .991568229456347;lower .921629512291253','research_v1/ld/native_coloc_all_signal_pairs.tsv','Native coloc5.2.3 all12 conditions','Known-region model-conditional shared association','Fixed priors/model;variant universe;no unique causal variant','Independent BF configuration/prior recomputation<=3.45e-15'),
('C07','rs2431108/rs77960 EUR r²=.9954806028','research_v1/ld/chr5_prior_vs_brain6_lead_LD.tsv','Actual ALT-dosage LD','Different lead IDs strongly tag the same region','No exact causal identity from LD','Public reference samples;not disease replication'),
('C08','Four clinical slots:2 pass/1 nonpass/1 unestimated','06_REPLICATION_RESULTS.tsv','Maximum native trait P;Bonferroni4','PHENOTYPE_SENSITIVITY; independence UNVERIFIED','Changed phenotype;filtering;cohort/control linkage','11 independent raw-JSON/arithmetic checks;no certified two-trait replication'),
('C09','Nominal5% illustrative null rejection4.08–52.16%','research_v1/statistics/calibration_simulations.tsv','Locked240000 truth-known draws','Source-specific uncertainty counterexample','Synthetic scale/LD;not observed empirical false-positive rate','12 unchanged-source fixtures;independent variance derivation'),
('C10','723/2199 derived r out of bounds;1105 missing','research_v1/statistics/empirical_block_diagnostics.tsv','All3304 inherited estimates audit','Historical covariance candidates have uncertified calibration','No corrected empirical SE/P;derived r not valid correlation claim','Direct numeric/determinant audit'),
('C11','16/72 conservative conditional contrasts','research_v1/global_contrasts/contrasts.tsv','Unknown covariance bound;Bonferroni72','Source-defined differences under marginal-normal assumptions','Retrospective;ascertainment/cohort/h²calibration;not mechanism','All72 independent stdlib tail/CI checks'),
('C12','0/60 all-variant molecular H4>=.8;max.61757469','research_v1/molecular/sensitivity_results.tsv','Five priors xthree QTLscale xfour inputs','No effector admitted from retained molecular contexts','Single-signal;source scale/strand;component-only','60 Python explicit-BF and60 source-Rfunction checks'),
('C13','0/60 nonpal molecular H4>=.8;max.45525975','research_v1/statistics/molecular_palindromic_diagnostic/sensitivity_results.tsv','Separate frozen exclusion stress grid','No robust high-H4 molecular result in reviewed grid','No biological absence inference;filteredrawlocal-only','Input hashes/pair harmonization local integration'),
('C14','22/25 protected windows have prior regional evidence','research_v1/NOVELTY_MASTER.tsv','Accessible primary supplements/build audit','Previously reported regions;threePD unresolved','Focused queries;some alleles/fulltext inaccessible','27 exact source/build checks;14 source-accounting checks'),
('C15','No candidate-C strong-QTL admission','research_v1/molecular/candidate_C_historical_QTL_admission.tsv','28 historical broader-window contexts','Molecular follow-up held;no causal gene','Gate failure not null;source/coverage restrictions','Copied exact eligibility rows;window containment'),
('C16','NO_GO original biological discovery','PUBLICATION_READINESS.md','Scientific contribution vsjournal criteria','Bounded investigator-review study;not submission-approved','Human review/permissions/independence/method advance outstanding','Automated checks distinct from human peer review')]
rows=[]
for key,result,path,method,wording,conditions,review in claims:
 p=A/path;assert p.is_file(),p
 rows.append({'claim_id':key,'result':result,'exact_input':path,'input_sha256':sha(p),'method':method,
    'allowed_manuscript_wording':wording,'failure_conditions_or_limits':conditions,'independent_evidence':review})
pd.DataFrame(rows).to_csv(A/'CLAIM_TO_EVIDENCE.tsv',sep='\t',index=False)
print(json.dumps({'current_source_rows':len(d),'components':len(components),'claims':len(rows),'candidate_C_QTL_contexts':len(c)}))
