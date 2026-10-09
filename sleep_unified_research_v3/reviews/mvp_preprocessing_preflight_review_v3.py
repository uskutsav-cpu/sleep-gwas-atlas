#!/usr/bin/env python3
"""Independent software preflight only; never materializes GWAS or estimates h2/rg."""
import ast
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parent
SSD = Path('/Volumes/Extreme SSD')
OLD = SSD/'Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'


def digest(path):
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    after = path.stat()
    return {'path': str(path), 'bytes': after.st_size, 'sha256': h.hexdigest(),
            'unchanged_size_mtime_during_hash': (before.st_size, before.st_mtime_ns) ==
                                              (after.st_size, after.st_mtime_ns)}


def main():
    acquisition_path = ROOT/'sleep_unified_research_v1/logs/mvp_insomnia_acquisition_receipt_v1.json'
    acquisition = json.loads(acquisition_path.read_text())
    paths = {
        'source': Path(acquisition['path']),
        'hm3_reference': SSD/'sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz',
        'hm3_alleles': OLD/'work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr/w_hm3.snplist',
        'chain': OLD/'work/track_b_completion/local_dependency_copies/ref/hg38ToHg19.over.chain.gz',
        'historical_47': ROOT/'discovery_extension/scripts/47_stream_replication_sources.py',
        'historical_streaming_helper': ROOT/'discovery_extension/scripts/streaming_io.py',
        'chain_helper': ROOT/'scripts/liftover_chain.py',
        'worker': PACKAGE/'scripts/19_materialize_mvp_hm3.py',
        'runner': PACKAGE/'scripts/20_run_mvp_preprocessing.py',
        'protocol': PACKAGE/'FROZEN_MVP_PREPROCESSING_PROTOCOL_v4.md',
        'python_binary': (OLD/'.ldsc-env/bin/python').resolve(),
        'acquisition_receipt': acquisition_path,
        'original_schema_receipt': ROOT/'sleep_unified_research_v1/logs/mvp_insomnia_source_schema_v1.json',
        'definition_contract': ROOT/'sleep_unified_research_v2/source_definition_review/provenance_v2_definition_contract.json',
        'protocol_amendment': ROOT/'sleep_unified_research_v2/manifests/protocol_amendment_v3_receipt.json',
        'independent_review_code': Path(__file__),
    }
    # The parent signaled that materialization started during review. Avoid a
    # concurrent large source read: verify its stat/header here and defer the
    # independent full source hash until the bounded worker has completed.
    hashes = {name: digest(path) for name, path in paths.items() if name != 'source'}
    hashes['source'] = {'path':str(paths['source']), 'bytes':paths['source'].stat().st_size,
                        'sha256_from_existing_acquisition_receipt':acquisition['actual_sha256'],
                        'independently_rehashed_this_preflight':False}
    expected = {
        'source': '46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835',
        'hm3_reference': 'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8',
        'hm3_alleles': 'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed',
        'chain': '14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1',
    }
    checks = []
    def check(name, actual, expected_value):
        checks.append({'name': name, 'pass': actual == expected_value,
                       'actual': actual, 'expected': expected_value})
    for name, value in expected.items():
        if name != 'source':
            check(name+'_pinned_actual_sha256', hashes[name]['sha256'], value)
    check('source_receipt_registered_sha256', acquisition['actual_sha256'], expected['source'])
    check('source_exact_bytes', hashes['source']['bytes'], 575504178)
    check('all_independently_hashed_files_stable', all(x['unchanged_size_mtime_during_hash'] for name,x in hashes.items() if name != 'source'), True)

    # Execute only the historical orientation function and its three constants.
    tree = ast.parse(paths['historical_47'].read_text())
    nodes = [node for node in tree.body if
             (isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and
              t.id in {'VALID', 'AMBIGUOUS', 'COMPLEMENT'} for t in node.targets)) or
             (isinstance(node, ast.FunctionDef) and node.name == 'orientation')]
    orientation_ns = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<orientation-function-only>', 'exec'), orientation_ns)
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    canonical = set(complement)
    ambiguous = [{'A','T'}, {'C','G'}]
    allele_cases = 0
    allele_failures = []
    for a1, a2 in itertools.permutations('ACGT', 2):
        if {a1,a2} in ambiguous:
            continue
        for effect, other in itertools.product(['A','C','G','T','N','AT',''], repeat=2):
            for strand in ['+', '-']:
                valid = effect in canonical and other in canonical and effect != other and {effect,other} not in ambiguous
                # Oracle treats the effect base and other base as identities;
                # strand conversion changes their names, never their roles.
                oracle = None
                if valid:
                    e = complement[effect] if strand == '-' else effect
                    o = complement[other] if strand == '-' else other
                    equivalence = [(e,o), (complement[e],complement[o])]
                    if (a1,a2) in equivalence:
                        oracle = 1
                    elif (a2,a1) in equivalence:
                        oracle = -1
                e = ''.join(complement.get(base,base) for base in effect) if strand == '-' else effect
                o = ''.join(complement.get(base,base) for base in other) if strand == '-' else other
                actual = orientation_ns['orientation'](e,o,a1,a2)
                allele_cases += 1
                if actual != oracle:
                    allele_failures.append([effect,other,a1,a2,strand,actual,oracle])
    check('allele_orientation_cases_no_mismatch', allele_failures, [])

    # Compare the imported point mapper to a separate exhaustive interval oracle.
    # These are small software fixtures, not biological or GWAS findings.
    spec = importlib.util.spec_from_file_location('reviewed_chain_helper', paths['chain_helper'])
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    intervals = [
        (10,15,2,100,1000,'+',1), (20,25,2,300,1000,'-',2),
        (30,35,2,400,1000,'+',3), (30,35,2,400,1000,'+',4),
        (40,45,2,500,1000,'+',5), (40,45,3,500,1000,'+',6),
        (50,55,None,600,1000,'+',7),
        (60,65,2,700,1000,'+',8), (60,65,2,295,1000,'-',9),
        (70,100,2,800,1000,'+',10), (80,81,2,900,1000,'+',11),
    ]
    index = helper.ChainIndex({1:list(intervals)})
    chain_failures = []
    for point in range(1,106):
        hits = set()
        for begin,end,chrom,start,size,strand,_ in intervals:
            if begin <= point-1 < end:
                displacement = point-1-begin
                # One-based endpoint formulas derived directly from UCSC.
                destination = start+displacement+1 if strand == '+' else size-start-displacement
                hits.add((chrom,destination,strand))
        if not hits:
            expected_point = ('unmapped',None)
        elif len(hits)>1:
            expected_point = ('ambiguous',None)
        elif next(iter(hits))[0] is None:
            expected_point = ('non_autosomal_target',None)
        else:
            expected_point = ('mapped',next(iter(hits)))
        actual_point = index.map_point('chr1',str(point))
        if actual_point != expected_point:
            chain_failures.append([point,actual_point,expected_point])
    check('105_chain_boundary_gap_overlap_points', chain_failures, [])
    for chromosome,position in [('chrX',1),(1,0),(1,-1),(1,'1.5'),(1,'bad')]:
        check('invalid_point_'+str(chromosome)+'_'+str(position), index.map_point(chromosome,position), ('invalid_source_coordinate',None))
    check('absent_autosome_unmapped', index.map_point(22,1), ('unmapped',None))

    # Fraction supplies a separate exact arithmetic expression for N_eff.
    exact_neff = Fraction(4*78566*329572,78566+329572)
    worker_neff = 4/(1/78566+1/329572)
    check('N_eff_exact_rational_agreement_within_one_float_ULP', abs(worker_neff-float(exact_neff)) <= math.ulp(float(exact_neff)), True)
    printed_neff = format(worker_neff,'.12g')
    check('N_eff_12_significant_digits', printed_neff, '253768.615047')
    ci_checks = []
    for beta,se in [(-.6,.1),(0,.2),(.4,.05)]:
        critical = 1.959963984540054
        odds = math.exp(beta)
        lower = math.exp(beta-critical*se)
        upper = math.exp(beta+critical*se)
        derived = math.log(odds)/((math.log(upper)-math.log(lower))/(2*critical))
        ci_checks.append(abs(derived-beta/se) <= 1e-12)
    check('CI_Z_log_scale_Wald_algebra_3_software_fixtures', all(ci_checks), True)

    # Check the reviewed pre-launch gates without running either worker or runner.
    worker = paths['worker'].read_text()
    runner = paths['runner'].read_text()
    tokens = {
        'exact_source_header': 'header!=expected_header' in worker,
        'unexpected_nonmissing_SE_fails': "elif supplied is None: raise RuntimeError('UNEXPECTED_NONMISSING_INVALID_SUPPLIED_SE')" in worker,
        'each_row_fixed_case_control_counts': '(408138,78566,329572)' in worker,
        'whole_file_19703815_row_gate': "counts['source_rows']!=19703815" in worker,
        'exclusive_row_reconciliation': "counts['source_rows']==sum(excluded.values())+counts['output_rows']" in worker,
        'first_eligible_duplicate_policy': worker.index("if rs in seen") > worker.index("if not math.isfinite(z)"),
        'final_GRCh37_MHC_inclusive': '25000000<=bp37<=34000000' in worker,
        'fresh_output_and_plan_files': "tmp.open('xb')" in worker and "path.open('x')" in runner and 'PRIOR_DESTINATION_PRESERVED' in runner,
        'monitor_partial_output': "Path(str(out)+'.partial')" in runner,
        'owned_worker_kill_escalation': 'proc.kill(); proc.wait()' in runner,
        'postexit_size_floor_time': all(x in runner for x in ['OUTPUT_SIZE_LIMIT_AFTER_EXIT','INTERNAL_EMERGENCY_FLOOR_AFTER_EXIT','TIME_LIMIT_AFTER_EXIT']),
        'input_hashes_pre_post': 'before==after' in runner,
        'frozen_plan_hash_bound_to_worker': "numerical['execution_plan_sha256']==frozen_plan_sha256" in runner and 'plan_unchanged' in runner,
        'no_h2_or_rg_admission': "'h2_or_rg_estimated':False" in runner and "'pair_tests_admitted':0" in runner,
    }
    for name,value in tokens.items():
        check('inspected_gate_'+name, value, True)
    with gzip.open(paths['source'],'rt',encoding='utf-8',newline='') as source:
        observed_header = source.readline().rstrip('\r\n').split('\t')
    schema = json.loads(paths['original_schema_receipt'].read_text())
    check('independent_current_source_header_matches_audited_header', observed_header, schema['source_header'])

    result = {
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'BOUNDED_PREFLIGHT_PASS' if all(item['pass'] for item in checks) else 'PREFLIGHT_FAILURE',
        'reviewer_role': 'separate independent numerical reproduction reviewer',
        'review_python': sys.version,
        'input_and_code_hashes': hashes,
        'hash_buffer_bytes':65536,
        'checks':checks,
        'check_count':len(checks),
        'allele_orientation_software_cases':allele_cases,
        'chain_geometry_software_cases':105,
        'N_eff_exact_fraction':str(exact_neff),
        'N_eff_exact_fraction_rounded_to_double':float(exact_neff),
        'N_eff_float_difference_from_exact_fraction':worker_neff-float(exact_neff),
        'N_eff_double':worker_neff,
        'N_eff_serialized_12g':printed_neff,
        'worker_or_runner_launched':False,
        'GWAS_materialization_performed':False,
        'h2_or_rg_estimated':False,
        'pair_tests_admitted':0,
        'official_chain_format_source':'https://genome.ucsc.edu/goldenPath/help/chain.html',
        'limits':[
            'Code inspection began before launch; parent launched the worker during this review after methods approval. Small software checks and sealing finish afterward without inspecting outcomes.',
            'Actual full source hash/stream, retained output, CRC and resource behavior await execution receipts and independent streaming review. Full source hashing deferred to avoid concurrent large reads.',
            'CI-derived Z is conditional on a 95% log-scale Wald interval; accession-specific uncertainty and effect semantics remain unconfirmed.',
            'Sampled resource thresholds are not hard instantaneous limits. Failure leaves logs and output preserved; an early worker failure may omit a full worker diagnostic receipt.',
            'No complete raw-chain reproduction, h2/rg, cohort disjointness, 396/1200-family completion or independent biological replication is certified.',
        ],
    }
    target = HERE/'mvp_preprocessing_preflight_review_v3.json'
    with target.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),
                      'allele_cases':allele_cases,'chain_cases':105,
                      'N_eff':worker_neff,'worker_sha256':hashes['worker']['sha256'],
                      'runner_sha256':hashes['runner']['sha256'],
                      'protocol_sha256':hashes['protocol']['sha256']}))
    if not all(item['pass'] for item in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
