#!/usr/bin/env python3
"""Independent, bounded, read-only arithmetic audit; no estimator or raw GWAS run.

Uses only Python stdlib. Writes new diagnostics only under this v4 directory.
Printed log precision bounds all calculations inherited from logs. The sample-
size convention sensitivity does not validate prevalence or native refit identity.
"""
import ast
import csv
import datetime
import gzip
import hashlib
import json
import math
import platform
import re
import statistics
import subprocess
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
FROZEN = REPO / 'sleep_unified_research_v1/sources/recovered'
ARCHIVE = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
MUNGED = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged')
TRAITS = ('ms', 'melanoma', 'ldl', 'hdl', 'triglycerides')
LIPIDS = ('ldl', 'hdl', 'triglycerides')
NUM = r'[-+\d.eE]+'
INPUTS = {}
CHECKS = []


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def bind(path):
    path = path.resolve()
    if str(path) not in INPUTS:
        INPUTS[str(path)] = dict(path=str(path), bytes=path.stat().st_size, sha256_before=sha(path))
    return path


def check(name, ok, evidence):
    CHECKS.append(dict(check=name, pass_=bool(ok), evidence=evidence))
    if not ok:
        raise AssertionError(name + ': ' + str(evidence))


def table(path):
    with bind(path).open() as f:
        return [dict(row, source_line=i) for i, row in enumerate(csv.DictReader(f, delimiter='\t'), 2)]


def content(path):
    return bind(path).read_text()


def value(pattern, text):
    m = re.search(pattern, text, re.M)
    if not m:
        raise ValueError('missing pattern: ' + pattern)
    return m


def summary(text):
    h = value(r'^Total (Observed|Liability) scale h2:\s*(' + NUM + r')\s*\((' + NUM + r')\)', text)
    i = value(r'^Intercept:\s*(' + NUM + r')\s*\((' + NUM + r')\)', text)
    result = dict(scale=h[1], h2=float(h[2]), h2_se=float(h[3]), h2_token=h[2], h2_se_token=h[3],
                  intercept=float(i[1]), intercept_se=float(i[2]), intercept_token=i[1], intercept_se_token=i[2])
    for label, pattern in [('mean_chi2', r'^Mean Chi\^2:\s*(' + NUM + r')'), ('lambda_gc', r'^Lambda GC:\s*(' + NUM + r')')]:
        result[label] = float(value(pattern, text)[1])
    m = re.search(r'^Ratio:\s*(' + NUM + r')', text, re.M)
    result['printed_ratio'] = float(m[1]) if m else None
    return result


def printed_half_unit(token):
    """Half of the last displayed decimal place, including scientific notation."""
    mantissa, *exponent = token.lower().split('e')
    decimals = len(mantissa.split('.')[1]) if '.' in mantissa else 0
    return .5 * 10 ** ((int(exponent[0]) if exponent else 0) - decimals)


def scan_n(trait, expected_hash):
    path = bind(MUNGED / (trait + '.sumstats.gz'))
    check('current_munged_hash_' + trait, INPUTS[str(path)]['sha256_before'] == expected_hash,
          'v1 native_input_hash_checks current-file hash; not a historical munged hash')
    rows = valid = 0
    nmin = math.inf
    nmax = -math.inf
    seen = set()
    total = 0.0
    correction = 0.0
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        check('munged_header_' + trait, header == ['SNP', 'A1', 'A2', 'Z', 'N'], header)
        for line in f:
            rows += 1
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 5:
                raise ValueError('malformed line ' + str(rows + 1))
            try:
                n, z = float(fields[4]), float(fields[3])
            except ValueError:
                continue
            if not (math.isfinite(n) and math.isfinite(z)):
                continue
            if n <= 0:
                raise ValueError('nonpositive N')
            valid += 1
            nmin, nmax = min(nmin, n), max(nmax, n)
            if len(seen) < 101:
                seen.add(n)
            # Compensated summation avoids misleading accumulation drift for constant N.
            y = n - correction
            t = total + y
            correction = (t - total) - y
            total = t
    return dict(trait=trait, path=str(path), source_sha256=INPUTS[str(path)]['sha256_before'],
                header='/'.join(header), full_template_rows=rows, finite_N_and_Z_rows=valid,
                missing_or_nonfinite_N_or_Z_rows=rows-valid, N_min=nmin, N_max=nmax,
                N_mean=total/valid, distinct_N_count_capped_at_101=len(seen),
                distinct_N_count_is_lower_bound=(len(seen) == 101), gzip_read_to_EOF=True,
                scope='processed source; no dense raw-source chain replay')


def conversion(p, k):
    if not (0 < p < 1 and 0 < k < 1):
        raise ValueError('invalid prevalence')
    t = statistics.NormalDist().inv_cdf(1-k)
    density = math.exp(-t*t/2) / math.sqrt(2*math.pi)
    return k*k*(1-k)*(1-k)/(p*(1-p)*density*density)


def write_table(name, rows):
    out = HERE / name
    if out.exists():
        raise FileExistsError('refusing to replace diagnostic: ' + str(out))
    with out.open('x') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    return dict(path=str(out), bytes=out.stat().st_size, sha256=sha(out), rows=len(rows))


def main():
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    bind(Path(__file__))
    panel_rows = table(FROZEN / 'config/analysis_panel.tsv')
    panel = {r['trait_id']: r for r in panel_rows}
    h2_rows = table(FROZEN / 'results/tables/h2_summary.tsv')
    h2 = {r['trait']: r for r in h2_rows}
    rg = table(FROZEN / 'results/tables/rg_matrix.tsv')
    hashes = table(REPO / 'sleep_unified_research_v1/tables/native_input_hash_checks.tsv')
    expected = {r['trait_id']: r['actual_sha256'] for r in hashes if r['kind'] == 'core_munged'}
    check('frozen_cardinality', len(panel_rows) == len(panel) == len(h2_rows) == 45 and len(rg) == 396,
          [len(panel_rows), len(panel), len(h2_rows), len(rg)])
    check('frozen_unique_pairs', len({(r['sleep_trait'], r['disease_trait']) for r in rg}) == 396,
          '12 sleep by 33 outcomes')
    code_path = REPO.parent / 'ldsc-code/ldscore/sumstats.py'
    code = content(code_path)
    archived_code = ARCHIVE / 'ldsc/ldscore/sumstats.py'
    check('pinned_control_flow_matches_archive', sha(bind(archived_code)) == sha(code_path), str(code_path))
    reg_code = content(REPO.parent / 'ldsc-code/ldscore/regressions.py')
    archived_reg = bind(ARCHIVE / 'ldsc/ldscore/regressions.py')
    check('pinned_regression_matches_archive', sha(archived_reg) == sha(REPO.parent / 'ldsc-code/ldscore/regressions.py'), str(archived_reg))
    for p in ['scripts/01_harmonize.py', 'scripts/03_h2_qc.sh', 'scripts/05_collate.py', 'scripts/25_genomicsem_covariance.R',
              'sleep_unified_research_v1/manifests/ldsc_native_code_comparison_v1.json']:
        bind(REPO / p)
    tree = ast.parse(code)
    funcs = {x.name: x for x in tree.body if isinstance(x, ast.FunctionDef)}
    split = funcs['_split_or_none']
    # Execute only this pure list-conversion function, never import LDSC.
    isolated = ast.Module(body=[split], type_ignores=[])
    namespace = {}
    exec(compile(isolated, str(code_path), 'exec'), namespace)
    converted = namespace['_split_or_none'](None, 34)
    fn = funcs['estimate_rg']
    branch = next(x for x in fn.body if isinstance(x, ast.If) and
                  'n_annot == 1' in ast.unparse(x.test) and 'args.two_step' in ast.unparse(x.test))
    prior_call = next(x for x in fn.body if isinstance(x, ast.Assign) and 'args.intercept_h2, args.intercept_gencov' in ast.unparse(x))
    active = eval(compile(ast.Expression(branch.test), str(code_path), 'eval'), {},
                  {'n_annot': 1, 'args': SimpleNamespace(two_step=None, intercept_h2=converted)})
    check('rg_default_branch_is_false', not active and converted == [None]*34 and prior_call.lineno < branch.lineno,
          dict(normalization_line=prior_call.lineno, branch_line=branch.lineno, branch_test=ast.unparse(branch.test)))
    hfn = funcs['estimate_h2']
    hbranch = next(x for x in ast.walk(hfn) if isinstance(x, ast.If) and
                   ast.unparse(x.test) == 'args.two_step is None and args.intercept_h2 is None')
    hactive = eval(compile(ast.Expression(hbranch.test), str(code_path), 'eval'), {},
                   {'args': SimpleNamespace(two_step=None, intercept_h2=None)})
    check('standalone_h2_default_branch_is_true', hactive, hbranch.lineno)
    control = dict(pinned_sumstats_sha256=sha(code_path), pinned_regressions_sha256=sha(REPO.parent/'ldsc-code/ldscore/regressions.py'),
                   normalization_line=prior_call.lineno, rg_branch_line=branch.lineno,
                   rg_branch_test=ast.unparse(branch.test), normalized_intercept_none_type=type(converted).__name__,
                   default_rg_automatic_two_step_active=active, standalone_h2_automatic_two_step_active=hactive,
                   h2_branch_line=hbranch.lineno, native_difference='UNESTIMATED; no estimator launched or code patched')

    n_rows = [scan_n(t, expected[t]) for t in TRAITS]
    n_by = {r['trait']: r for r in n_rows}
    standalone = {}
    for trait in TRAITS:
        path = ARCHIVE/'results/logs'/('h2_'+trait+'.log')
        text = content(path)
        s = summary(text)
        s.update(trait=trait, path=str(path), two_step_message='Using two-step estimator with cutoff at 30.' in text,
                 input_rows=int(value(r'^Read summary statistics for (\d+) SNPs', text)[1]),
                 merged_rows=int(value(r'^After merging with regression SNP LD, (\d+) SNPs remain', text)[1]),
                 h2_log_line=next(i for i,l in enumerate(text.splitlines(),1) if l.startswith('Total ')),
                 intercept_log_line=next(i for i,l in enumerate(text.splitlines(),1) if l.startswith('Intercept:')))
        check('h2_finite_source_count_' + trait, s['input_rows'] == n_by[trait]['finite_N_and_Z_rows'], s['input_rows'])
        check('standalone_two_step_' + trait, s['two_step_message'], str(path))
        standalone[trait] = s
    cache = {}
    lipid_rows = []
    pair_rows = []
    log_modes = []
    for r in rg:
        if r['disease_trait'] not in TRAITS:
            continue
        logpath = ARCHIVE / r['input_log']
        if str(logpath) not in cache:
            text = content(logpath)
            segments = list(re.finditer(r'^Computing rg for phenotype \d+/\d+', text, re.M))
            parsed = {}
            for i,m in enumerate(segments):
                end = segments[i+1].start() if i+1 < len(segments) else len(text)
                segment = text[m.start():end]
                other = value(r'Reading summary statistics from .*?/([^/\s]+)\.sumstats\.gz', segment)[1]
                if other not in TRAITS:
                    continue
                offset = segment.index('Heritability of phenotype ')
                if segment[offset:].startswith('Heritability of phenotype 1\n'):
                    offset = segment.index('Heritability of phenotype ', offset+1)
                hsection = segment[offset:segment.index('Genetic Covariance', offset)]
                s = summary(hsection)
                s['log_segment_line'] = text[:m.start()].count('\n')+1
                s['h2_log_line'] = text[:m.start()+offset].count('\n')+1+next(j for j,l in enumerate(hsection.splitlines()) if l.startswith('Total '))
                s['intercept_log_line'] = text[:m.start()+offset].count('\n')+1+next(j for j,l in enumerate(hsection.splitlines()) if l.startswith('Intercept:'))
                s['valid_allele_SNPs'] = int(value(r'^(\d+) SNPs with valid alleles\.', segment)[1])
                s['source_full_template_rows'] = int(value(r'^Read summary statistics for (\d+) SNPs', segment)[1])
                parsed[other] = s
            cache[str(logpath)] = parsed
            log_modes.append(dict(path=str(logpath), two_step_message='Using two-step estimator' in text,
                                  two_step_cli='--two-step' in text, explicit_intercept_cli='--intercept-h2' in text))
        s = cache[str(logpath)][r['disease_trait']]
        trait = r['disease_trait']
        for field, key in [('h2_obs','h2'),('h2_obs_se','h2_se'),('h2_int','intercept'),('h2_int_se','intercept_se')]:
            # frozen collator uses 4 significant digits; authoritative log may print more.
            check('log_table_printed_' + str(r['source_line']) + '_' + field,
                  float(r[field]) == float(format(s[key], '.4g')), dict(table=r[field], log=s[key]))
        check('pair_full_template_count_' + str(r['source_line']), s['source_full_template_rows'] == n_by[trait]['full_template_rows'], s['source_full_template_rows'])
        detail = dict(table_line=r['source_line'], sleep_trait=r['sleep_trait'], disease_trait=trait,
                      input_log=str(logpath), log_segment_line=s['log_segment_line'], h2_log_line=s['h2_log_line'],
                      intercept_log_line=s['intercept_log_line'], rg=float(r['rg']), rg_se=float(r['se']), p=float(r['p']),
                      frozen_FDR_396=float(r['fdr']), frozen_analysis_tier=r['analysis_tier'],
                      historical_outcome_verdict=r['disease_h2_verdict'], pair_h2_N_input_scale=s['h2'],
                      pair_h2_SE_N_input_scale=s['h2_se'], pair_h2_int=s['intercept'], pair_h2_int_SE=s['intercept_se'],
                      pair_intercept_over_1_2=s['intercept']>1.2, pair_valid_allele_SNPs=s['valid_allele_SNPs'])
        pair_rows.append(detail)
        if trait in LIPIDS:
            detail = dict(detail)
            detail.update(standalone_h2=standalone[trait]['h2'], standalone_h2_SE=standalone[trait]['h2_se'],
                          standalone_intercept=standalone[trait]['intercept'], standalone_intercept_SE=standalone[trait]['intercept_se'],
                          standalone_SNPs=standalone[trait]['merged_rows'],
                          pair_minus_standalone_intercept=s['intercept']-standalone[trait]['intercept'],
                          pair_intercept_distance_to_1_2=s['intercept']-1.2,
                          descriptive_pair_int_95CI_low=s['intercept']-1.959963984540054*s['intercept_se'],
                          descriptive_pair_int_95CI_high=s['intercept']+1.959963984540054*s['intercept_se'],
                          pair_mean_chi2=s['mean_chi2'], pair_lambda_gc=s['lambda_gc'],
                          diagnostic_inflation_ratio=(s['intercept']-1)/(s['mean_chi2']-1),
                          source_N_min=n_by[trait]['N_min'], source_N_max=n_by[trait]['N_max'],
                          native_two_step_effect='UNESTIMATED', current_inference='HOLD_UNQUALIFIED_LIPID_CLAIM',
                          limit='shared sources and changed SNP sets; no SE for intercept difference without covariance; descriptive CI does not override gate')
            lipid_rows.append(detail)
    check('target_pair_counts', len(lipid_rows) == 36 and len(pair_rows) == 60,
          dict(lipids=len(lipid_rows), targeted_all=len(pair_rows)))
    check('lipid_all_pair_intercepts_above_existing_boundary', all(r['pair_intercept_over_1_2'] for r in lipid_rows), '36/36')
    check('default_pair_logs_no_two_step', all(not r['two_step_message'] and not r['two_step_cli'] and not r['explicit_intercept_cli'] for r in log_modes), log_modes)

    liability = []
    pair_liability = []
    for trait in ('ms','melanoma'):
        p = panel[trait]
        s = standalone[trait]
        text = content(Path(s['path']))
        p_logged = float(value(r'--samp-prev (' + NUM + r')', text)[1])
        k_logged = float(value(r'--pop-prev (' + NUM + r')', text)[1])
        nc, nt = int(p['ncase']), int(p['ncontrol'])
        total = nc + nt
        sample = nc/total
        neff = 4/(1/nc+1/nt)
        input_n = n_by[trait]['N_min']
        check('constant_effective_N_' + trait, input_n == n_by[trait]['N_max'] and abs(input_n-neff) < .00051,
              dict(input_N=input_n, exact_Neff=neff))
        check('logged_prevalence_' + trait, p_logged == round(sample,6) and k_logged == float(p['pop_prev']), [p_logged,k_logged])
        c_original = conversion(p_logged,k_logged)
        c_balanced = conversion(.5,k_logged)
        factor = c_balanced/c_original
        check('balanced_factor_cancels_K_' + trait, math.isclose(factor,4*p_logged*(1-p_logged),rel_tol=2e-15), factor)
        original_fit = s['h2']/c_original
        original_se = s['h2_se']/c_original
        balanced_h2 = original_fit*c_balanced
        balanced_se = original_se*c_balanced
        total_route_h2 = original_fit*input_n/total*conversion(sample,k_logged)
        total_route_se = original_se*input_n/total*conversion(sample,k_logged)
        check('fixed_fit_Z_invariance_' + trait, math.isclose(s['h2']/s['h2_se'],balanced_h2/balanced_se,rel_tol=2e-15), s['h2']/s['h2_se'])
        row = dict(trait=trait, panel_line=p['source_line'], h2_summary_line=h2[trait]['source_line'],
                   original_log=s['path'], original_h2_line=s['h2_log_line'], cases=nc, controls=nt, actual_total_N=total,
                   actual_case_fraction=sample, logged_case_fraction=p_logged, source_N=input_n, exact_Neff=neff,
                   frozen_population_K=k_logged, original_reported_liability_h2=s['h2'], original_reported_liability_SE=s['h2_se'],
                   original_conversion_factor=c_original, balanced_conversion_factor=c_balanced,
                   original_fit_h2_on_Neff_input_scale=original_fit, original_fit_SE_on_Neff_input_scale=original_se,
                   N_total_scale_h2_fixed_fit=original_fit*input_n/total, N_total_scale_SE_fixed_fit=original_se*input_n/total,
                   balanced_same_fit_liability_h2=balanced_h2, balanced_same_fit_liability_SE=balanced_se,
                   total_N_actual_P_route_fixed_fit_liability_h2=total_route_h2,
                   total_N_actual_P_route_fixed_fit_liability_SE=total_route_se,
                   route_relative_difference_due_to_N_rounding=(total_route_h2-balanced_h2)/balanced_h2,
                   original_to_balanced_multiplier=1/factor, rescaling_factor=factor,
                   unchanged_Z=s['h2']/s['h2_se'], original_intercept=s['intercept'],
                   historical_QC_verdict=h2[trait]['verdict'], diagnostic_Z_ge_4=s['h2']/s['h2_se']>=4,
                   original_log_Neff_h2_proxy=neff*s['h2'], balanced_same_fit_Neff_h2_proxy=neff*balanced_h2,
                   inherited_h2_half_print_unit=printed_half_unit(s['h2_token']),
                   balanced_h2_half_print_unit=printed_half_unit(s['h2_token'])*factor,
                   inherited_SE_half_print_unit=printed_half_unit(s['h2_se_token']),
                   balanced_SE_half_print_unit=printed_half_unit(s['h2_se_token'])*factor,
                   provenance='rounded original standalone log and direct current processed N; no native rerun',
                   scientific_status='SCALE_INCONSISTENCY_DIAGNOSED; PREVALENCE_AND_ASCERTAINMENT_NOT_RESOLVED')
        liability.append(row)
        for r in pair_rows:
            if r['disease_trait'] != trait:
                continue
            pair_liability.append(dict(table_line=r['table_line'],sleep_trait=r['sleep_trait'],disease_trait=trait,
                                       pair_h2_input_Neff_scale=r['pair_h2_N_input_scale'],pair_h2_input_Neff_SE=r['pair_h2_SE_N_input_scale'],
                                       balanced_K_unchanged_liability_h2=r['pair_h2_N_input_scale']*c_balanced,
                                       balanced_K_unchanged_liability_SE=r['pair_h2_SE_N_input_scale']*c_balanced,
                                       unchanged_h2_Z=r['pair_h2_N_input_scale']/r['pair_h2_SE_N_input_scale'],
                                       historical_rg=r['rg'],historical_rg_SE=r['rg_se'],historical_p=r['p'],historical_FDR_396=r['frozen_FDR_396'],
                                       scale_presentation_changes_rg=False, scale_presentation_changes_intercept=False,
                                       native_N_recoding_effect='UNESTIMATED', historical_tier=r['frozen_analysis_tier'],
                                       status='ARITHMETIC_SENSITIVITY_ONLY; not new outcomes or new QC admission'))

    sem_rows = table(ARCHIVE/'results/tables/ldsc_covariance_pairs.tsv')
    sem_inclusion = table(ARCHIVE/'results/tables/genomicsem_trait_inclusion.tsv')
    sem_diag = []
    for r in sem_rows:
        if r['trait_1'] != r['trait_2'] or r['trait_1'] not in TRAITS:
            continue
        t = r['trait_1']
        fraction = int(panel[t]['ncase'])/(int(panel[t]['ncase'])+int(panel[t]['ncontrol'])) if t in ('ms','melanoma') else None
        factor = 4*fraction*(1-fraction) if fraction is not None else 1
        inclusion = next(x for x in sem_inclusion if x['trait_id']==t)
        sem_diag.append(dict(trait=t, covariance_table_line=r['source_line'], inclusion_table_line=inclusion['source_line'],
                             original_diagonal=float(r['genetic_covariance']), original_diagonal_SE=float(r['genetic_covariance_se']),
                             original_intercept=float(r['cross_trait_intercept']), original_include=inclusion['include_genomic_sem'],
                             scale_factor_if_original_actual_P_input=factor,
                             diagnostic_rescaled_diagonal=float(r['genetic_covariance'])*factor,
                             diagnostic_rescaled_SE=float(r['genetic_covariance_se'])*factor,
                             limit='conditional diagonal arithmetic only; not a corrected S/V export or calibrated correlation uncertainty'))
    primary = [r for r in rg if r['analysis_tier']=='PRIMARY_PHASE1' and float(r['fdr'])<.05]
    all_fdr = [r for r in rg if float(r['fdr'])<.05]
    summary_counts = dict(core_panel_rows=len(panel_rows),core_h2_table_rows=len(h2_rows),core_rg_rows=len(rg),
                          primary_positive_under_frozen_396=sum(1 for r in primary), all_positive_under_frozen_396=len(all_fdr),
                          lipid_rows=len(lipid_rows),lipid_historical_fdr_positive=sum(r['frozen_FDR_396']<.05 for r in lipid_rows),
                          lipid_positive_by_trait={t:sum(r['frozen_FDR_396']<.05 for r in lipid_rows if r['disease_trait']==t) for t in LIPIDS},
                          lipid_intercept_ranges={t:[min(r['pair_h2_int'] for r in lipid_rows if r['disease_trait']==t),max(r['pair_h2_int'] for r in lipid_rows if r['disease_trait']==t)] for t in LIPIDS},
                          legacy_primary_positive_without_20_flagged_lipid_rows=sum(r['disease_trait'] not in LIPIDS for r in primary),
                          direct_processed_source_scans=5,liability_traits=2,liability_pair_rows=len(pair_liability),
                          original_GenomicSEM_elements=len(sem_rows),original_GenomicSEM_traits=len(sem_inclusion))
    check('historical_frozen_family_counts', summary_counts['primary_positive_under_frozen_396']==153 and summary_counts['all_positive_under_frozen_396']==161,
          summary_counts)
    artifacts = [write_table('lipid_intercept_diagnostics_v4.tsv',lipid_rows),
                 write_table('processed_N_coding_diagnostics_v4.tsv',n_rows),
                 write_table('liability_input_sensitivity_v4.tsv',liability),
                 write_table('liability_pair_scale_sensitivity_v4.tsv',pair_liability),
                 write_table('legacy_GenomicSEM_diagonal_sensitivity_v4.tsv',sem_diag)]
    bind(HERE/'official_lipid_readme_fetch_receipt_v4.json')
    for acc in ('GCST90239658','GCST90239652','GCST90239664'):
        bind(HERE/(acc+'_official_README.txt'))
    for r in INPUTS.values():
        r['sha256_after']=sha(Path(r['path']))
        check('input_unchanged_'+r['path'],r['sha256_before']==r['sha256_after'],r['path'])
    receipt = dict(schema='independent_statistical_arithmetic_v4.1',started_utc=started,
                   completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=platform.python_version(),
                   git_branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),
                   git_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
                   status='DIAGNOSTIC_ARITHMETIC_PASS_SCIENTIFIC_QC_UNRESOLVED',summary=summary_counts,
                   control_flow=control,standalone_original_logs=list(standalone.values()),pairwise_log_modes=log_modes,
                   inputs=list(INPUTS.values()),outputs=artifacts,checks=CHECKS,
                   no_native_estimator_called=True,no_dense_raw_source_read=True,no_QC_threshold_changed=True,
                   original_testing_family_unchanged=True,prior_packages_unchanged=True,
                   limitations=['Direct N audit covers only MS, melanoma, LDL, HDL, and triglycerides; other binary N modes may be variant-specific.',
                                'Current processed hashes match current ledger, not an unavailable historical processed-file hash.',
                                'All original h2/log arithmetic inherits printed precision; output digits do not recover original full precision.',
                                'P=.5 is an N_eff convention sensitivity, not an assertion of biological sample prevalence.',
                                'No alternative population prevalence selected or invented; legacy K is held fixed.',
                                'No native explicit two-step comparison or N recoding rerun; effect on native estimates unestimated.',
                                'Changing liability presentation rescales h2 and its SE after fit and leaves Z, intercept, and rg unchanged.',
                                'An N recoding rerun can change estimation through iterative weights/clipping/filtering; algebra is conditional on a fixed fit.',
                                'No calibrated GenomicSEM covariance matrix or shared-estimator heterogeneity covariance established.'])
    out=HERE/'qc_arithmetic_receipt_v4.json'
    with out.open('x') as f:
        json.dump(receipt,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(dict(status=receipt['status'],checks=len(CHECKS),summary=summary_counts,receipt=str(out),receipt_sha256=sha(out)),indent=2))


if __name__=='__main__':
    main()
