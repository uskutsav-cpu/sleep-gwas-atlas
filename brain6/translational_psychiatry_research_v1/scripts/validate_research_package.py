"""Read-only cross-artifact validation; does not certify biological inference."""
from pathlib import Path
import hashlib,json,re,math
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
A=ROOT/'brain6/translational_psychiatry_research_v1';R=A/'research_v1'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
checks=[]
def check(name,condition,detail):
    checks.append({'check':name,'status':'PASS' if bool(condition) else 'FAIL','detail':detail})
def read(p):return pd.read_csv(p,sep='\t')
claims=read(A/'CLAIM_TO_EVIDENCE.tsv')
check('Complete16claim_source_bindings',len(claims)==16 and all(sha(A/row.exact_input)==row.input_sha256 for row in claims.itertuples()),'Exact result artifact hashes, methods, limits and independent evidence fields')
g=read(A/'MANUSCRIPT/source_tables/BRAIN6_FINAL_GLOBAL.tsv')
check('Inherited72_35_original_family',len(g)==72 and g.significance_under_original_396_family.sum()==35,'No subset FDR recomputation')
lava=read(A/'MANUSCRIPT/source_tables/BRAIN6_FINAL_LAVA.tsv')
check('Immutable_LAVA_failure',len(lava)==17465 and lava.status.eq('TESTED').sum()==13745 and lava.status.eq('NOT_RUN').sum()==3720 and lava.family_decision.eq('FAILED_QC_NOT_PROMOTED').all(),'NOT_RUN3720 exceeds locked873; not biological null')
c=read(A/'MANUSCRIPT/source_tables/BRAIN6_FINAL_CANDIDATES.tsv')
check('Candidate_count',len(c)==25,'Pair-specific protected screening candidates')
fit=read(R/'ld/native_finemap_fit_summary.tsv')
check('Eight_actual_native_fits',len(fit)==8 and fit.converged.all() and fit.n_snps.eq(2185).all() and fit.n_credible_sets.eq(1).all(),'Same-discovery exploratory reference-qualified fits')
co=read(R/'ld/native_coloc_all_signal_pairs.tsv')
primary=co[co.model.eq('PRIMARY_L10_MEDIAN_NEFF')]
h4=lambda prior:float(primary.loc[(primary.p12-prior).abs()<1e-12,'PP.H4.abf'].iloc[0])
check('Native12posterior_conditions',len(co)==12 and ((co[[f'PP.H{i}.abf' for i in range(5)]].sum(axis=1)-1).abs()<1e-10).all(),'Every signal-pair/prior result retained')
check('Primary_H4_values',math.isclose(h4(1e-5),.991568229456347,abs_tol=1e-12) and math.isclose(h4(1e-6),.921629512291253,abs_tol=1e-12),'Default and prespecified lower prior; no causal identity inference')
clinical=read(A/'06_REPLICATION_RESULTS.tsv')
check('Four_slot_clinical_outcomes',len(clinical)==4 and clinical.decision.eq('BOTH_ASSOCIATIONS_AT_FIXED_VARIANT').sum()==2 and clinical.decision.eq('NOT_ESTIMATED').sum()==1 and clinical.decision.eq('DID_NOT_MEET_FAMILY_THRESHOLD').sum()==1,'Two pass;one nonpass;one unestimated under Bonferroni4')
check('Clinical_scope_not_promoted',clinical.classification.eq('PHENOTYPE_SENSITIVITY').all() and clinical.cohort_independence.eq('UNVERIFIED').all() and clinical.same_causal_variant.eq('NOT_TESTED').all(),'No independent two-trait/shared-causal-variant claim')
cal=read(R/'statistics/calibration_simulations.tsv')
check('Locked_truth_known_calibration',len(cal)==24 and cal.replicates.sum()==240000 and cal.classification.eq('TRUTH_KNOWN_SIMULATION').all(),'Not empirical false-positive calibration')
check('Recorded_rejection_range',math.isclose(cal['p_alpha0.05_rate'].min(),.0408) and math.isclose(cal['p_alpha0.05_rate'].max(),.5216),'Truth-null nominal5% scenario fractions')
emp=read(R/'statistics/empirical_block_diagnostics.tsv')
check('All_historical_covariance_estimates',len(emp)==3304 and emp.derived_correlation_out_of_bounds.sum()==723 and emp['corr'].isna().sum()==1105 and emp.source_family_significant.sum()==3,'Historical arithmetic, uncertified model calibration')
summ=json.loads((R/'global_contrasts/summary.json').read_text())
check('Conditional72contrast_family',summ['planned_contrasts']==72 and summ['bounded_differences']==16 and len(read(R/'global_contrasts/covariance_sensitivity.tsv'))==648,'Retrospective same-source normal-model diagnostic')
mol=read(R/'molecular/sensitivity_results.tsv');pal=read(R/'statistics/molecular_palindromic_diagnostic/sensitivity_results.tsv')
check('All120molecular_conditions_negative',len(mol)==len(pal)==60 and max(mol.pp_h4.max(),pal.pp_h4.max())<.8,'No effector admitted; not proof of biological absence')
nov=read(R/'NOVELTY_MASTER.tsv');protected=nov[nov.candidate_locus_id.astype(str).str.startswith('placo_')]
check('Novelty_protected25_known22',len(protected)==25 and protected.classification.eq('KNOWN_SAME_PAIR_REGION').sum()==22,'Three PD unresolved rather than novel')
gate=read(R/'molecular/candidate_C_historical_QTL_admission.tsv')
check('CandidateC_QTL_admission_hold',len(gate)==28 and not (pd.to_numeric(gate.min_qtl_p_shared,errors='coerce')<=5e-8).any(),'Broader historical window contains new native C; no admitted strong QTL')
source=read(A/'SOURCE_AND_ACCESS_LEDGER.tsv')
check('Fresh_sources_and_source_specific_terms',len(source)==27 and source.provenance_scope.str.contains('CURRENT').sum()>=12 and source.license.astype(str).str.contains('forbids public reposting').any(),'Historical rows and fresh recovery/reference/FinnGen feasibility explicit; four molecular resources retained')
tex=(A/'MANUSCRIPT/brain6_reassessment.tex').read_text()
abstract=tex.split('\\section*{Abstract}')[1].split('\\textbf{Keywords:}')[0]
strip=lambda s:re.sub(r'\\[A-Za-z]+(?:\[[^\]]*\])?', '',re.sub(r'\\(?:cite|href)\{[^}]+\}', '',s)).replace('{',' ').replace('}',' ')
abswords=len(strip(abstract).split());mpwords=len((A/'MANUSCRIPT/MP_ABSTRACT.txt').read_text().split())
check('Abstracts_and_qualified_conclusions',abswords<=350 and 150<=mpwords<=250 and 'previously reported chr5' in abstract and 'preclude a new causal or mechanistic discovery claim' in abstract,'GM words='+str(abswords)+';MP words='+str(mpwords))
body=tex.split('\\section*{Background}')[1].split('\\section*{Abbreviations}')[0]
check('Conservative_article_length',len(strip(body).split())<=3500,'Main body estimated words='+str(len(strip(body).split())))
check('Verified_reference_keys',len(re.findall(r'\\bibitem\{',tex))==15 and len(read(A/'MANUSCRIPT/REFERENCES_VERIFIED.tsv'))==15,'15 references;primary metadata/access limitations retained')
check('Public_figure_outputs',all((A/'MANUSCRIPT/figures'/f'figure{n}_{name}.{ext}').is_file() for n,name in [(1,'global_map'),(2,'native_chr5'),(3,'covariance_diagnostic')] for ext in ['pdf','svg','png']),'Three vector figures with source tables and script/hash provenance')
bad=[str(p.relative_to(A)) for p in A.rglob('*') if p.is_file() and (p.suffix.lower()=='.rds' or p.name.endswith('_signed_LD.npz') or ('nonpal' in p.name and p.suffix=='.gz'))]
check('No_restricted_native_or_filtered_artifacts_in_public_area',not bad,str(bad))
allowed={'COMPLETE_VALIDATED','COMPLETE_EXPLORATORY','FAILED_QC','NOT_RUN','BLOCKED_DATA','BLOCKED_ACCESS','BLOCKED_SOFTWARE','NEEDS_REVIEW','NOT_APPLICABLE'}
inventory=read(A/'COMPONENT_INVENTORY.tsv')
check('Component_status_schema',len(inventory)==28 and set(inventory.status)<=allowed,'Scientific evidence/software/human holds classified separately')
result={'classification':'DERIVED_ARTIFACT_AND_MANUSCRIPT_VALIDATION_NOT_SCIENTIFIC_APPROVAL',
    'checks':checks,'passed':sum(c['status']=='PASS' for c in checks),'failed':sum(c['status']=='FAIL' for c in checks)}
print(json.dumps(result,indent=2))
if result['failed']:raise SystemExit(1)
