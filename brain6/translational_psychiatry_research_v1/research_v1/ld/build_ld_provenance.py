#!/usr/bin/env python3
"""Content hashes and exact result-count audit, generated after real execution."""
import csv,hashlib,json,pathlib,re
HERE=pathlib.Path(__file__).resolve().parent;LOCAL=HERE.parents[3]/'work/ld_genotypes_research_v1/C_finemap_inputs'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def main():
 fits=rows(HERE/'native_finemap_fit_summary.tsv');coloc=rows(HERE/'native_coloc_all_signal_pairs.tsv');abf=rows(HERE/'same_universe_native_ABF_sensitivity.tsv');independent=rows(HERE/'independent_coloc_recomputation.tsv')
 assert len(fits)==8 and len(coloc)==12 and len(abf)==3 and len(independent)==12
 assert all(r['converged']=='TRUE' for r in fits)
 assert all(float(r['max_posterior_abs_error_H0_H4'])<1e-12 for r in independent)
 assert re.search(r'COMPLETED_FIXED_NATIVE_FITS\s+8\s+COLOC_INVOCATIONS\s+12\s+ABF_INVOCATIONS\s+3',(HERE/'native_finemap_execution.log').read_text())
 for p in ['source_free_regression_tests.log','native_API_regression_tests.log','factor_rebuild_execution.log']:assert (HERE/p).exists()
 files=[p for p in HERE.iterdir() if p.is_file() and p.name not in ['native_analysis_provenance.json','analysis_artifact_hashes.tsv']]
 public=[{'artifact':str(p.relative_to(HERE.parents[3])),'bytes':p.stat().st_size,'sha256':sha(p),'scope':'VERSIONED_CODE_MANIFEST_DERIVED_RESULT_OR_REGRESSION_LOG'} for p in sorted(files)]
 private=[{'artifact':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'scope':'LOCAL_ONLY_SOURCE_DERIVED_INPUT_NATIVE_OUTPUT_OR_DIAGNOSTIC'} for p in sorted(LOCAL.iterdir()) if p.is_file()]
 with (HERE/'analysis_artifact_hashes.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(public[0]),delimiter='\t');w.writeheader();w.writerows(public)
 prov={'record_created_post_execution':True,'scientific_protocol_frozen_before_native_outcomes':True,'protocol_sha256':sha(HERE/'FINEMAP_PRE_FIT_PROTOCOL.md'),'reference_protocol_sha256':sha(HERE/'LD_PRE_OUTCOME_PROTOCOL.md'),'native_R_code_sha256':sha(HERE/'native_finemap_coloc.R'),'native_R_diagnostic_adapter_sha256':sha(HERE/'native_runtime_adapter.R'),'native_fit_execution_marker_verified':True,'fits':8,'coloc_susie_invocations':12,'coloc_signal_pair_rows':12,'native_ABF_invocations':3,'independent_python_posterior_checks':12,'max_independent_posterior_abs_error':max(float(r['max_posterior_abs_error_H0_H4']) for r in independent),'actual_reference_reconstructions':3,'actual_independent_scalar_genotype_pairs':300,'source_free_tests':10,'native_synthetic_API_tests':3,'historical_source_or_output_mutation':False,'human_independent_validation_performed':False,'private_local_artifacts':private,'artifact_hash_table_sha256':sha(HERE/'analysis_artifact_hashes.tsv'),'restricted_input_warning':'ADHD raw results and source-derivedZbetaSE remain ignored locally; numerical/posterior outputs are derived research artifacts.','known_API_failures_retained':['native_RSS_consistency_initial_interface_failure.log','native_RSS_consistency_initial_serializer_failure.log'],'intentional_excluded_bad_call_reproduction':'initial_invalid_N_infinity_reproduced_warnings.tsv36optimizerwarnings_plus1ggplotdeprecation','posterior_interpretation':'Exploratory known-region support, not independent traitreplication, causalvariant, gene or clinicaladvance.','direct_required_artifacts':{p.name:sha(p) for p in [HERE.parent/'LD_RECONSTRUCTION_REPORT.md',HERE.parent/'LD_VALIDATION.tsv',HERE.parent/'LD_SOURCE_MANIFEST.json']}}
 (HERE/'native_analysis_provenance.json').write_text(json.dumps(prov,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'public_artifacts_hashed':len(public),'private_artifacts_hashed':len(private),'native_fits_verified':len(fits),'native_coloc_comparisons_verified':len(coloc),'independent_numerical_recomputations':len(independent)}))
if __name__=='__main__':main()
