"""Write research evidence pointers and environment/checksum records after execution."""
import datetime, json, platform, sys
import numpy as np
import pandas as pd
import scipy
import matplotlib
import pytest
from calibration import HERE, sha, write_tsv, verify_freeze

def main():
    verify_freeze()
    environment={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'executable':sys.executable,'python':sys.version,'platform':platform.platform(),'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'pytest':pytest.__version__},'execution_limits':{'processes':1,'batch_size':128,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},'native_scope':'Unchanged upstream calLocalCov function on synthetic LD API only; no genotype/GWAS pipeline.'}
    (HERE/'statistics_environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    claims=[
        ('STAT_SOURCE_IDENTITY','historical_source_identity.tsv','6/6 historical Python file hashes match','SOURCE_IDENTITY_VALIDATED','Exact implementation identity; no genotype provenance proof'),
        ('STAT_NATIVE_CONCORDANCE','upstream_concordance.tsv','12 unchanged upstream numerical fixtures','SYNTHETIC_ROUTINE_CONCORDANCE','Synthetic LD API; no full native empirical integration'),
        ('STAT_CALIBRATION_COUNTEREXAMPLE','calibration_simulations.tsv','24 scenarios;240000 replicates;0 invalid;5% rejection4.08–52.16%','IMPLEMENTATION_SPECIFIC_MODEL_COUNTEREXAMPLE','No empirical false-positive estimate or variance replacement'),
        ('STAT_NORMAL_TAIL_LIMIT','gaussian_product_tail_diagnostics.tsv','Effective moment variance modes1.35–112.10;excess kurtosis0.054–4.456','POST_RESULT_DETERMINISTIC_DIAGNOSTIC','Correct variance does not imply exact normal family tail'),
        ('STAT_EMPIRICAL_BOUNDS','empirical_pair_diagnostics.tsv','3304 estimates;723/2199 numeric correlations out of bounds;1105 missing','RETROSPECTIVE_SOURCE_DIAGNOSTIC','No bounded r claim;no joint h2/covariance determinant test'),
        ('STAT_EMPIRICAL_FRAGILITY','significant_block_fragility.tsv','3 reported FWER positives;SE multipliers1.2371/1.0538/1.4011 lose threshold','HISTORICAL_REPORTED_POSITIVE_UNCERTIFIED_CALIBRATION','Frozen denominator8465;no calibrated corrected empirical P'),
        ('STAT_SCALAR_N_INVARIANCE','scalar_n_invariance.tsv','4 fixed-Z scalar N rescalings;Wald P and r invariant','ALGEBRAIC_IDENTITY','Unchanged P cannot validate per-variant N/model semantics'),
        ('STAT_REGRESSION_TESTS','regression_tests.log','9 meaningful source-free numerical/inference regression tests passed','SOFTWARE_NUMERICAL_VALIDATION','No empirical GWAS or independent cohort validation'),
        ('STAT_PARENT_GLOBAL_REVIEW','parent_global_contrast_checks.tsv','All72contrasts independently recomputed;16conditional differences;48unique marginalrows','RETROSPECTIVE_SOURCE_PARAMETER_CONTRAST_REVIEW','Conditional on calibrated jointnormal margins;not16independent biologicaldiscoveries'),
        ('STAT_PARENT_MOLECULAR_REVIEW','parent_molecular_configuration_checks.tsv','All60conditions independently recomputed using explicit distinct-variant configuration sum;maxerror3.45e-15','ACTUAL_INHERITED_SLICE_NUMERICAL_REVIEW','Same4ADHDcomponentcontexts;no insomnia molecular result;no independentreplication'),
        ('STAT_REFERENCE_R_REVIEW','reference_R_posterior_checks.tsv','60retainedcoloc5.2.3referencefunction calculations;maxposteriorerror8.39e-15','REFERENCE_R_FUNCTION_EXECUTION','Fullcolocpackage orchestration notinstalledor rerun;sdYassumptionspersist'),
        ('STAT_DENSE_VARIANCE_REVIEW','parent_independent_variance_checks.tsv','2stackedGaussiandensequadraticformchecks;maxrelativeerror6.54e-16;parentMCwithin0.540MCSE','POST_RESULT_INDEPENDENT_NUMERICAL_REVIEW','Corroborates idealmodel identity;notempiricalGWAScalibration'),
        ('STAT_PALINDROMIC_STRESS','molecular_palindromic_diagnostic/sensitivity_results.tsv','60newfrozenconditions;allsourcegatespass;maxH4.455260;zeroH4>=.8','POST_RESULT_HARMONIZATION_SENSITIVITY','Palindromicdoesnotproveambiguousidentity;SVmodel/sdY/onecomponentconstraintsremain'),
        ('STAT_PARENT_REVIEW_TESTS','parent_review_regression_tests.log','7additional source-free inference/output regression tests passed','SOFTWARE_NUMERICAL_VALIDATION','Separatefromoriginal9tests;zeroindependentcohortreplication'),
        ('STAT_LOCAL_INTEGRATION_TESTS','parent_review_local_integration_tests.log','2actuallocalintegrationtests passed:restrictedfilteredmolecularslicesandnativefit/BF/LDbyteidentity','ACTUAL_LOCAL_INPUT_IDENTITY_INTEGRATION','Missinginputs fail when requested;not independentcohortreplication or fullfits rerun'),
        ('STAT_NATIVE_FINE_MAPPING_REVIEW','native_fine_mapping_review_summary.json','8savednativefits/8credible sets and12signalpairpriorconditions independentlyreviewed;primaryH4.991568/lowerprior.921630','SAME_SOURCE_NATIVE_MODEL_NUMERICAL_REVIEW','Knownchr5region;oneCSpertrait;binarymedian-N/external503referenceapproximation;no newcausalgene or independentreplication')
    ]
    rows=[{'claim_id':i,'method':'SUPERGNOVA_source_specific_audit','hypothesis_family':'Historical8465 preserved;new24 truth-known scenarios prospectively frozen','evidence_file':p,'evidence_sha256':sha(HERE/p),'result':r,'classification':c,'limitations':l,'independent_two_trait_replications':0} for i,p,r,c,l in claims]
    write_tsv('STATISTICAL_EVIDENCE_LEDGER.tsv',rows)
    files=[p for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='STATISTICS_CHECKSUMS.sha256']
    (HERE/'STATISTICS_CHECKSUMS.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(HERE)}\n' for p in files))

if __name__=='__main__': main()
