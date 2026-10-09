#!/usr/bin/env python3
"""Small receipt/capture command binding audit; never executes a native fit."""
import datetime
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
spec = importlib.util.spec_from_file_location('independent_reader', Path(__file__).with_name('independent_core_subset_checker_v4_3.py'))
U = importlib.util.module_from_spec(spec);spec.loader.exec_module(U)


def main():
    U.text(Path(__file__))
    plan = U.source_json(P/'manifests/ssd_native_execution_plan_v4_3.json')
    support = Path(plan['ssd_support_package'])
    native = support/'native/core_reproduction_v1'
    python = next(s for s in plan['dependencies_sha256'] if s.endswith('/.ldsc-env/bin/python'))
    ref = str(Path(next(s for s in plan['dependencies_sha256'] if s.endswith('/1.l2.ldscore.gz'))).parent)+'/'
    ldsc_dir = str(ROOT.parent/'ldsc-code')
    proofs = []
    for job in (j for j in plan['jobs'] if j['stage'] == 'core'):
        prefix = native/job['job_id']
        r = U.source_json(Path(str(prefix)+'.execution_receipt.json'))
        c = U.source_json(Path(str(prefix)+'.full_precision.json'))
        expected = [python,'-u',str(support/'scripts/native_ldsc_capture.py'),'--ldsc-dir',ldsc_dir,
                    '--'+job['kind'],','.join(job['inputs']),'--ref-ld-chr',ref,'--w-ld-chr',ref,
                    '--print-delete-vals','--out',str(prefix)]+job['options']
        args = c['arguments']
        opts = dict(zip(job['options'][::2],job['options'][1::2]))
        checks = {'exact_command':r['command'] == expected, 'capture_ldsc_dir':c['ldsc_dir'] == ldsc_dir,
                  'capture_input_argument':args[job['kind']] == ','.join(job['inputs']),
                  'capture_output_prefix':args['out'] == str(prefix),
                  'capture_ref_prefix':args['ref_ld_chr'] == ref and args['w_ld_chr'] == ref,
                  'capture_print_deletes':args['print_delete_vals'] is True and args['n_blocks'] == 200,
                  'samp_prev_bound_to_job':args['samp_prev'] == opts.get('--samp-prev'),
                  'pop_prev_bound_to_job':args['pop_prev'] == opts.get('--pop-prev'),
                  'allele_checks_enabled':args['no_check_alleles'] is False,
                  'free_intercepts':args['no_intercept'] is False and args['intercept_h2'] is None and args['intercept_gencov'] is None,
                  'reference_M_5_50_default':args['not_M_5_50'] is False,
                  'capture_libraries':c['libraries'] == {'numpy':'1.21.5','pandas':'1.3.3','scipy':'1.7.3'},
                  'capture_python':c['python'].startswith('3.9.23'),
                  'original_outputs_unmodified_claim':r['original_outputs_modified'] is False,
                  'outputs_bounded_to_job_namespace':all(Path(s).parent == native and Path(s).name.startswith(job['job_id']) and not Path(s).name.startswith('._') for s in r['all_output_sha256'])}
        proofs.append({'job_id':job['job_id'],'checks':checks,'all_pass':all(checks.values())})
    unchanged = all(U.sha(Path(s)) == m['sha256'] for s,m in U.SEEN.items())
    out = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'57core receipt/capture command,prevalence,reference,environment bindings only; no native execution',
           'all_pass':len(proofs) == 57 and all(r['all_pass'] for r in proofs) and unchanged,
           'jobs':proofs,'inputs':U.SEEN,'consumed_hashes_unchanged':unchanged,'original_no_mutation_claim_is_not_full_filesystem_audit':True}
    dest = P/'reviews/independent_core_command_binding_receipt_v4.json'
    assert not dest.exists();dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'job_count':len(proofs),'all_pass':out['all_pass']},indent=2))
    assert out['all_pass']


if __name__ == '__main__':
    main()
