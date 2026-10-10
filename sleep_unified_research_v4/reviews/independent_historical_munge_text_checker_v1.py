#!/usr/bin/env python3
"""Independent, synthetic-only qualification of the pre-existing header-text munger.

No real GWAS/reference bodies are opened. All actual CLI inputs are invented below.
Original modules are read only; subprocesses use -B and owned SSD TMPDIR.
"""
import ast
import bz2
import csv
import difflib
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time

P = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OWN = SSD / 'tmp/independent_historical_munge_text_fixture_v1_2'
OUT = P / 'reviews/independent_historical_munge_text_receipt_v1.json'
PROOF = P / 'source_provenance/historical_munge_header_text_compatibility_v1.json'
PIN = P.parent.parent / 'ldsc-code'
T0 = time.monotonic()
BINDINGS = {}
CONTROLS = []
RUNS = []
PEAK = 0

def sha(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(65536), b''):
            h.update(b)
    return h.hexdigest()

def bind(path, expected=None):
    path = Path(path)
    got = sha(path)
    if expected is not None:
        assert got == expected, (str(path), got, expected)
    BINDINGS[str(path)] = {'sha256': got, 'bytes': path.stat().st_size}
    return got

def check(name, condition, detail=None):
    CONTROLS.append({'name': name, 'pass': bool(condition), 'detail': detail})
    assert condition, (name, detail)

def budget():
    assert time.monotonic() - T0 < 480
    assert shutil.disk_usage(P).free >= 128 * 2**20
    assert shutil.disk_usage(SSD).free >= 2**30
    if OWN.exists():
        assert sum(p.stat().st_size for p in OWN.rglob('*') if p.is_file()) < 8 * 2**20

def write(path, data):
    with path.open('xb') as f:
        f.write(data)

def input_versions(name, text):
    data = text.encode('ascii')
    plain = OWN / (name + '.txt')
    gz = OWN / (name + '.txt.gz')
    bz = OWN / (name + '.txt.bz2')
    write(plain, data)
    write(gz, gzip.compress(data, mtime=0))
    write(bz, bz2.compress(data))
    return [plain, gz, bz]

def launch(name, command):
    global PEAK
    budget()
    env = dict(os.environ)
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', TMPDIR=str(OWN / 'tmp'))
    log = OWN / (name + '.stdout.log')
    start = time.monotonic()
    peak = 0
    with log.open('xb') as f:
        proc = subprocess.Popen(command, stdout=f, stderr=subprocess.STDOUT, env=env,
                                cwd=OWN, start_new_session=True)
        try:
            while proc.poll() is None:
                budget()
                rss = subprocess.run(['ps', '-o', 'rss=', '-p', str(proc.pid)],
                                     capture_output=True, text=True).stdout.strip()
                peak = max(peak, int(rss or '0') * 1024)
                assert peak < 256 * 2**20, ('RSS', peak)
                assert time.monotonic() - start < 30, ('deadline', name)
                time.sleep(.025)
            rc = proc.wait()
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
    PEAK = max(PEAK, peak)
    RUNS.append({'name': name, 'command': command, 'returncode': rc,
                 'stdout_sha256': sha(log), 'elapsed_seconds': time.monotonic() - start,
                 'peak_observed_worker_RSS_bytes': peak})
    return rc, log.read_text()

def cli(name, exe, source, extra):
    prefix = OWN / name
    rc, log = launch(name, [str(PYTHON), '-B', str(exe), '--sumstats', str(source),
                            '--out', str(prefix)] + extra)
    content = None
    rows = None
    target = Path(str(prefix) + '.sumstats.gz')
    if target.exists():
        # Invented file only; reaching EOF checks gzip CRC/trailer through stdlib.
        with gzip.open(target, 'rb') as f:
            content = f.read(65536)
            assert not f.read(1)
        rows = list(csv.DictReader(io.StringIO(content.decode('ascii')), delimiter='\t'))
    return rc, log, content, rows

def z_from_p(p):
    # Independently invert two-sided Gaussian probability, no scipy/stock helper.
    if p == 1:
        return 0.0
    lo, hi = 0.0, 40.0
    for _ in range(120):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def main():
    global PYTHON
    assert not OUT.exists()
    assert not OWN.exists(), 'Owned synthetic fixture target immutable; choose a new version explicitly.'
    assert shutil.disk_usage(P).free >= 256 * 2**20
    OWN.mkdir(parents=True)
    (OWN / 'tmp').mkdir()
    j = json.loads(PROOF.read_text())
    bind(PROOF, '8e9388c92d6801c299a1d61b7d59fa91da789fd2c2422f4f4983e721dbcdb41a')
    cand = Path(j['candidate_executable'])
    old = Path(j['pinned_unmodified_munger'])
    bind(old, j['pinned_unmodified_munger_sha256'])
    bind(cand, j['candidate_executable_sha256'])
    meta = Path(j['historical_metadata_path'])
    bind(meta, j['historical_metadata_sha256'])
    mj = json.loads(meta.read_text())
    check('historical_metadata_commit_dirty_status_and_candidate_hash',
          mj['commit'] == j['pinned_commit'] and mj['preexisting_dirty_status'] == 'M munge_sumstats.py'
          and any(x['sha256'] == sha(cand) and x['path'].endswith('/munge_sumstats.py')
                  for x in mj['source_files']))
    for path, receipt in j['source_and_copy_identity'].items():
        bind(path, receipt['sha256'])
        assert Path(path).stat().st_size == receipt['bytes']
    check('all_11_historical_source_and_copy_identities', len(j['source_and_copy_identity']) == 22)
    for path, expected in j['unchanged_pinned_imported_module_sha256'].items():
        bind(path, expected)
        rel = str(Path(path).relative_to(PIN))
        blob = subprocess.check_output(['git', '-C', str(PIN), 'show', j['pinned_commit'] + ':' + rel])
        check('pinned_commit_blob_' + rel, hashlib.sha256(blob).hexdigest() == expected)
    check('pinned_munger_commit_blob', hashlib.sha256(subprocess.check_output(
          ['git', '-C', str(PIN), 'show', j['pinned_commit'] + ':munge_sumstats.py'])).hexdigest() == sha(old))
    oldtext, newtext = old.read_text(), cand.read_text()
    diff = ''.join(difflib.unified_diff(oldtext.splitlines(True), newtext.splitlines(True),
                 fromfile='pinned_unmodified', tofile='previously_existing_historical_executable'))
    check('exact_declared_text_diff_only', diff == j['exact_diff'], diff)
    def strip_header(text):
        tree = ast.parse(text)
        found = [x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == 'read_header']
        assert len(found) == 1
        tree.body.remove(found[0])
        return ast.dump(tree, include_attributes=False)
    check('whole_scientific_AST_identical_except_read_header', strip_header(oldtext) == strip_header(newtext))
    plan = P / 'manifests/ssd_native_execution_plan_v4_3.json'
    bind(plan)
    pj = json.loads(plan.read_text())
    PYTHON = Path(next(k for k in pj['dependencies_sha256'] if '/.ldsc-env/bin/python' in k))
    bind(PYTHON, pj['dependencies_sha256'][str(PYTHON)])
    code = ('import sys,json,numpy,pandas,scipy,locale;print(json.dumps({'
            '"python":sys.version.split()[0],"numpy":numpy.__version__,'
            '"pandas":pandas.__version__,"scipy":scipy.__version__,'
            '"preferred_encoding":locale.getpreferredencoding(False)}))')
    rc, versionlog = launch('runtime_versions', [str(PYTHON), '-B', '-c', code])
    versions = json.loads(versionlog)
    check('pinned_runtime_versions', rc == 0 and {k: versions[k] for k in ['python','numpy','pandas','scipy']} ==
          {'python':'3.9.23','numpy':'1.21.5','pandas':'1.3.3','scipy':'1.7.3'}, versions)
    failed = SSD / 'core_pipeline/original_small_input_replay_v1'
    for rel in ['core_original_small_execution_receipt_v1.json','receipts_v4/ms__munge.worker.json',
                'logs_v4/ms__munge.worker.stdout.log']:
        bind(failed / rel)
    worker = json.loads((failed / 'receipts_v4/ms__munge.worker.json').read_text())
    failurelog = (failed / 'logs_v4/ms__munge.worker.stdout.log').read_text()
    check('preserved_production_failure_is_same_pinned_header_typeerror',
          worker['command'][2] == str(old) and worker['returncode'] == 1
          and worker['owned_cleanup_verified'] and "TypeError: a bytes-like object is required, not 'str'" in failurelog
          and 'read_header' in failurelog)
    text = 'rsid a1 a2 n p-value beta INFO eaf\r\n'
    data = [
      ['rpos','A','C',100,.05,.04,.9,.2], ['rneg','C','A',101,.01,-.04,1,.3],
      ['rone','G','T',150,1,-.02,1,.4], ['rsmall','a','c',200,1e-8,.06,1,.4],
      ['rzero','T','G',110,.5,0,1,.2], ['dup','A','G',140,.02,-.01,1,.4],
      ['dup','A','G',140,.02,.01,1,.4], ['rbadmerge','G','A',140,.1,.01,1,.4],
      ['rinfo','A','C',140,.1,.01,.899,.4], ['rinfoNA','A','C',140,.1,.01,'NA',.4],
      ['rmaf','A','C',140,.1,.01,1,.01], ['rmafhi','A','C',140,.1,.01,1,.99],
      ['rpzero','A','C',140,0,.01,1,.4], ['rpbad','A','C',140,1.1,.01,1,.4],
      ['rat','A','T',140,.1,.01,1,.4], ['rcg','C','G',140,.1,.01,1,.4],
      ['rindel','A','AC',140,.1,.01,1,.4], ['rnlow','A','C',99,.1,.01,1,.4],
      ['rbetaNA','A','C',140,.1,'NA',1,.4], ['rpNA','A','C',140,'NA',.01,1,.4],
      ['rnNA','A','C','NA',.1,.01,1,.4], ['.','A','C',140,.1,.01,1,.4],
      ['notref','A','C',140,.1,.01,1,.4]]
    text += ''.join(' '.join(map(str, row)) + '\r\n' for row in data)
    sources = input_versions('filters_alias_CRLF', text)
    alleles = OWN / 'invented_alleles.txt'
    order = ['rsmall','rpos','rneg','rone','rzero','dup','rbadmerge','absent']
    write(alleles, b'SNP A1 A2\nrsmall T G\nrpos A C\nrneg A C\nrone G T\nrzero T G\ndup A G\nrbadmerge A C\nabsent A C\n')
    flags = ['--merge-alleles', str(alleles), '--n-min', '100', '--chunksize', '4']
    candidate = []
    for fmt, source in zip(['plain','gzip','bz2'], sources):
        run = cli('candidate_' + fmt, cand, source, flags)
        check('candidate_' + fmt + '_CLI_success', run[0] == 0, run[1][-800:])
        candidate.append(run)
    runold = cli('pinned_plain', old, sources[0], flags)
    check('pinned_plain_CLI_success', runold[0] == 0, runold[1][-800:])
    check('plain_gzip_bz2_and_pinned_plain_identical_decompressed_output',
          all(x[2] == runold[2] for x in candidate), hashlib.sha256(runold[2]).hexdigest())
    for fmt, source in zip(['gzip','bz2'], sources[1:]):
        rc, log, content, _ = cli('pinned_' + fmt + '_failure', old, source, flags)
        check('preserved_pinned_' + fmt + '_bytes_str_failure', rc != 0 and content is None
              and "TypeError: a bytes-like object is required, not 'str'" in log)
    rows = candidate[0][3]
    check('merge_reference_order_and_missing_rows', [r['SNP'] for r in rows] == order)
    expected = {'rpos':(100,.05,.04,'A','C'), 'rneg':(101,.01,-.04,'C','A'),
                'rone':(150,1,-.02,'G','T'), 'rsmall':(200,1e-8,.06,'A','C'),
                'rzero':(110,.5,0,'T','G'), 'dup':(140,.02,-.01,'A','G')}
    for row in rows:
        if row['SNP'] not in expected:
            check('missing_merge_row_' + row['SNP'], all(v == '' for k,v in row.items() if k != 'SNP'))
            continue
        n,p,b,a1,a2 = expected[row['SNP']]
        z = z_from_p(p) * (-1 if b < 0 else 1)
        check('independent_N_sign_PZ_alleles_precision_' + row['SNP'],
              row['N'] == '%.3f' % n and row['Z'] == '%.3f' % z and row['A1'] == a1 and row['A2'] == a2,
              {'observed':row,'expected_Z':'%.3f'%z})
    # No merge, so all filter and duplicate dispositions can be observed directly.
    no_merge = cli('candidate_filter_accounting', cand, sources[1], ['--n-min','100','--chunksize','4'])
    check('filter_dispositions_and_duplicate_first', no_merge[0] == 0 and
          [r['SNP'] for r in no_merge[3]] == ['rpos','rneg','rone','rsmall','rzero','dup','rbadmerge','rmafhi','notref'],
          [r['SNP'] for r in no_merge[3]] if no_merge[3] else no_merge[1][-600:])
    check('stock_upper_literal_point99_behavior_preserved', 'rmafhi' in [r['SNP'] for r in no_merge[3]],
          'Stock min(f,1-f)>0.01 retains literal .99 due binary roundoff; no scientific threshold correction made.')
    check('no_invented_INFO_or_frequency_in_default_output', set(rows[0]) == {'SNP','A1','A2','N','Z'})
    keep = cli('candidate_keep_maf', cand, sources[1], ['--n-min','100','--keep-maf','--chunksize','4'])
    check('keep_maf_preserves_frequency_only_when_explicit', keep[0] == 0 and
          all('FRQ' in r and 'INFO' not in r for r in keep[3]))
    # Separate input designs exercise stock case/control and constant-N branches.
    case = input_versions('case_control', 'SNP A1 A2 P BETA N_CASES N_CONTROLS\n'
                          'cc1 A C .05 .02 50 150\ncc2 A C .01 -.02 25 75\ncc3 A C .2 0 30 170\n')
    caseflags = ['--n-min','1','--chunksize','2']
    cr = cli('candidate_case_control_gzip', cand, case[1], caseflags)
    co = cli('pinned_case_control_plain', old, case[0], caseflags)
    check('case_control_same_as_pinned_plain', cr[0] == 0 and co[0] == 0 and cr[2] == co[2])
    # Max-total-N rows cc1/cc3 have mean case fraction (.25+.15)/2=.20.
    check('stock_case_control_N_formula', [r['N'] for r in cr[3]] == ['250.000','125.000','150.000'],
          'N=(cases+controls)*(cases/total)/mean(case fraction among max-total-N rows); not 4*Ncas*Ncon/(Ncas+Ncon).')
    const = input_versions('OR_constant_N', 'SNP A1 A2 P OR\no1 A C .05 1.02\no2 G T .01 .98\no3 A G 1 1\n')
    cnflags = ['--N','200','--chunksize','2']
    cn = cli('candidate_OR_constant_N_bz2', cand, const[2], cnflags)
    cp = cli('pinned_OR_constant_N_plain', old, const[0], cnflags)
    check('OR_null1_fixedN_same_as_pinned_plain', cn[0] == 0 and cp[0] == 0 and cn[2] == cp[2]
          and [r['N'] for r in cn[3]] == ['200.000']*3
          and [r['Z'] for r in cn[3]] == ['%.3f'%z_from_p(.05),'%.3f'%-z_from_p(.01),'0.000'])
    default = input_versions('default_N', 'SNP A1 A2 N P BETA\nn1 A C 100 .05 .02\n'
                             'n2 A C 100 .01 -.02\nn3 A C 50 .2 0\n')
    nd = cli('candidate_default_N_gzip', cand, default[1], [])
    check('stock_default_N_90th_percentile_over_1point5', nd[0] == 0 and
          [r['SNP'] for r in nd[3]] == ['n1','n2'] and '66.66666666666667' in nd[1])
    for badname, badtext, error in [
      ('missing_P', 'SNP A1 A2 N BETA\nx A C 100 .01\n', 'Could not find P column.')]:
        b = input_versions(badname, badtext)
        ca = cli('candidate_' + badname, cand, b[1], [])
        po = cli('pinned_' + badname, old, b[0], [])
        check('unchanged_parser_rejection_' + badname, ca[0] != 0 and po[0] != 0
              and error in ca[1] and error in po[1])
    # Pin an existing parser limitation rather than asserting a nonexistent guard.
    duplicate = input_versions('duplicate_header', 'SNP A1 A2 N P P BETA\nx A C 100 .05 .05 .01\n')
    ca = cli('candidate_duplicate_header', cand, duplicate[1], [])
    po = cli('pinned_duplicate_header', old, duplicate[0], [])
    check('unchanged_legacy_duplicate_P_header_acceptance', ca[0] == 0 and po[0] == 0 and ca[2] == po[2],
          'Pinned indentation returns during first header-column loop iteration before later duplicate-P check; future exact canonical input-schema validator remains required.')
    correction = P / 'reviews/independent_historical_munge_text_initial_control_correction_v1.json'
    bind(correction)
    initial = json.loads(correction.read_text())
    bind(P / 'reviews/independent_historical_munge_text_initial_control_v1.py', initial['initial_checker_sha256'])
    for path, expected in initial['initial_fixture_file_sha256'].items():
        bind(path, expected)
    for path, receipt in BINDINGS.items():
        assert sha(path) == receipt['sha256'], ('mutation', path)
    budget()
    artifacts = {str(p): {'sha256':sha(p),'bytes':p.stat().st_size}
                 for p in OWN.rglob('*') if p.is_file()}
    result = {'schema':'independent_historical_munge_text_qualification_v1',
              'status':'QUALIFIED_PRELAUNCH_METHODS_PASS', 'controls':CONTROLS,
              'passed_control_count':len(CONTROLS), 'fixture_subprocesses':RUNS,
              'consumed_metadata_and_code_before_after_sha256':BINDINGS,
              'synthetic_fixture_artifacts':artifacts, 'runtime_versions':versions,
              'scientific_AST_identical_except_read_header':True,
              'native_fit_implementation_modified':False,
              'production_GWAS_body_reads':0,'production_workers':0,'estimator_calls':0,
              'historical_per_trait_binary_attestation':False,'scientific_full_replay_certified':False,
              'execution_admission':False,
              'qualification':'ASCII compressed header decoding only, prior dirty historical code at pinned commit. '
                              'Scientific parser/filter/N/P-to-Z/allele/missingness/3-decimal code unchanged. '
                              'Actual future production ordered content/QC agreement remains required.',
              'resource_bounds':{'worker_RSS_max_bytes':256*2**20,'worker_deadline_seconds':30,
                                 'stage_deadline_seconds':480,'owned_output_cap_bytes':8*2**20,
                                 'internal_launch_floor_bytes':256*2**20,'internal_emergency_floor_bytes':128*2**20},
              'elapsed_seconds':time.monotonic()-T0,'peak_observed_worker_RSS_bytes':PEAK,
              'reviewer_maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    with OUT.open('x') as f:
        json.dump(result,f,indent=2,sort_keys=True); f.write('\n')
    print(json.dumps({'status':result['status'],'controls':len(CONTROLS),'subprocesses':len(RUNS),
                      'receipt_sha256':sha(OUT),'peak_RSS':PEAK,'elapsed':result['elapsed_seconds']}))

if __name__ == '__main__':
    main()
