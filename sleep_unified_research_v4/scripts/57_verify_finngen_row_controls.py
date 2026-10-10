#!/usr/bin/env python3
"""Small result-free numerical/identity controls; never reads a GWAS body."""
import importlib.util
import json
from pathlib import Path
import sys

from canonical_calibration_common_v4_3 import sha, utc, write_new


def main():
    scripts = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('feasibility_preprocessor', scripts/'56_preprocess_finngen_insomnia_feasibility.py')
    pre = importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)
    root = scripts.parents[1]
    code = root/'scripts/liftover_chain.py'
    spec = importlib.util.spec_from_file_location('original_chain_controls', code)
    lift = importlib.util.module_from_spec(spec); spec.loader.exec_module(lift)
    chain = lift.ChainIndex({1:[(100, 110, 1, 200, 1000, '+', 1),
                               (120, 130, 1, 300, 1000, '-', 2)]})
    hm3 = {'rs1':('A','C',1,201), 'rs2':('T','G',1,700), 'rs3':('A','C',6,25_000_000),
           'rs4':('A','C',6,34_000_000), 'rs5':('A','C',6,34_000_001)}
    base = dict(zip(pre.HEADER, ['1','101','C','A','rs1','GENE','0.0455002638963584','1.34','0.2','0.1','0.3','0.3','0.3']))
    checks = []
    def check(label, actual, expected):
        checks.append(dict(label=label, actual=actual, expected=expected, passed=actual==expected))
        if actual != expected:
            raise AssertionError(label)
    def qualify(changes=None, lookup=None, mapper=None, seen=None):
        return pre.qualify_row({**base, **(changes or {})}, lookup or hm3, mapper or chain, seen or set())
    reason, row = qualify()
    check('forward_exact_coordinate_retained', reason, 'retained')
    check('forward_supplied_beta_se_Z', row['Z'], 2.0)
    check('effect_swap_signed_Z', qualify({'alt':'C','ref':'A'})[1]['Z'], -2.0)
    check('strand_complement_orientation', qualify({'alt':'T','ref':'G'})[1]['Z'], 2.0)
    reason, row = qualify({'pos':'121','rsids':'rs2'})
    check('reverse_chain_coordinate_retained', reason, 'retained')
    check('reverse_chain_orientation', row['Z'], 2.0)
    check('reverse_chain_flag', row['reverse_chain'], True)
    check('chain_1_based_positive', chain.map_point('1','101'), ('mapped',(1,201,'+')))
    check('chain_1_based_negative', chain.map_point('1','121'), ('mapped',(1,700,'-')))
    check('exclusive_interval_end', chain.map_point('1','111'), ('unmapped',None))
    check('invalid_source_coordinate', qualify({'pos':'0'})[0], 'liftover_invalid_source_coordinate')
    check('wrong_reference_coordinate', qualify({'pos':'102'})[0], 'mapped_coordinate_rsid_reference_mismatch')
    check('ambiguous_multiple_hm3_ids', qualify({'rsids':'rs1,rs2'})[0], 'multiple_hm3_rsids')
    check('duplicate_alias_same_identity', qualify({'rsids':'rs1;rs1'})[0], 'retained')
    check('strict_rsid_token_boundary', qualify({'rsids':'rs1x'})[0], 'not_nonambiguous_hm3')
    check('unknown_rsid', qualify({'rsids':'rs99'})[0], 'not_nonambiguous_hm3')
    check('retained_duplicate_rejected', qualify(seen={'rs1'})[0], 'duplicate_retained_hm3_rsid')
    check('ambiguous_alleles_rejected', qualify({'alt':'A','ref':'T'})[0], 'allele_mismatch_or_ambiguous')
    check('indel_rejected', qualify({'alt':'AT'})[0], 'allele_mismatch_or_ambiguous')
    check('maf_boundary_rejected', qualify({'af_alt':'0.01'})[0], 'maf_at_or_below_0_01_or_invalid')
    check('nonfinite_frequency_rejected', qualify({'af_alt':'nan'})[0], 'maf_at_or_below_0_01_or_invalid')
    check('zero_supplied_SE_rejected', qualify({'sebeta':'0'})[0], 'invalid_beta_or_supplied_se')
    check('negative_supplied_SE_rejected', qualify({'sebeta':'-1'})[0], 'invalid_beta_or_supplied_se')
    check('nonfinite_beta_rejected', qualify({'beta':'nan'})[0], 'invalid_beta_or_supplied_se')
    check('nonfinite_Z_rejected', qualify({'beta':'1e300','sebeta':'1e-300'})[0], 'nonfinite_beta_over_se')
    check('missing_P_does_not_filter', qualify({'pval':'NA'})[0], 'retained')
    check('inconsistent_P_does_not_filter', qualify({'pval':'1e-200'})[0], 'retained')
    check('literal_zero_P_does_not_filter', qualify({'pval':'0'})[0], 'retained')
    class MHCMapper:
        def __init__(self, p): self.p=p
        def map_point(self, c, p): return 'mapped',(6,self.p,'+')
    for rsid, position in [('rs3',25_000_000),('rs4',34_000_000),('rs5',34_000_001)]:
        check('MHC_coordinate_'+str(position), qualify({'rsids':rsid}, mapper=MHCMapper(position))[0],
              'extended_MHC_GRCh37' if position<=34_000_000 else 'retained')
    ambiguous = lift.ChainIndex({1:[(100,110,1,200,1000,'+',1),(100,110,1,201,1000,'+',2)]})
    check('ambiguous_chain_rejected', qualify(mapper=ambiguous)[0], 'liftover_ambiguous')
    receipt = dict(schema='result_free_finngen_row_controls_v1', completed_utc=utc(),
                   dependencies_sha256={str(p):sha(p) for p in [Path(__file__),scripts/'56_preprocess_finngen_insomnia_feasibility.py',code]},
                   checks=checks, all_pass=all(r['passed'] for r in checks),
                   GWAS_body_reads=0, estimator_calls=0, reference_full_body_reads=0,
                   fixture_scope='synthetic identity and supplied-statistic algebra controls only; no biological results or sampling calibration')
    target = root/'sleep_unified_research_v4/logs/finngen_row_controls_v4.json'
    write_new(target, receipt)
    print(json.dumps(dict(check_count=len(checks), all_pass=receipt['all_pass'], receipt=str(target))))


if __name__ == '__main__':
    main()
