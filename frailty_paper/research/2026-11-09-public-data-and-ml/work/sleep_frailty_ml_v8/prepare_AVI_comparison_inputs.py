"""Permissible composite AVI inputs for a future matched experiment, not labels."""
import sqlite3, collections, concurrent.futures, threading, time
from v8_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
import pysam

URL = 'https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
INDEX_MD5 = '911ee04b6e8fba1882465221b5a14f2c'


def main():
    assert_previous_intact()
    bound = json.loads((V7/'outputs/AVI_collection_receipt.json').read_text())
    for name, digest in bound['outputs'].items():
        assert sha(V7/name) == digest
    index = V7/'sources/AVI_SNVS.tsv.gz.tbi'
    assert hashlib.md5(index.read_bytes()).hexdigest() == INDEX_MD5
    with np.load(DATA/'measured_training_arrays.npz') as z:
        nb = int(z['bank_index'].max())+1
    training = np.load(V7/'avi_bank.npy')
    assert training.shape == (nb, 3)
    # Keep exact original index/allele alignment. The new dataset explicitly
    # retained V7 bank indices; no genotype or source effect is imputed here.
    np.save(DATA/'avi_training_bank.npy', training)
    wanted = {r['variant_GRCh38'] for r in csv_rows(DATA/'reserved_matoba_physical_design.csv') if r['source_baseline_allelic_coefficient_available'] == 'True'}
    old = sqlite3.connect((V7/'avi_exact_scores.sqlite').as_uri()+'?mode=ro', uri=True)
    db = sqlite3.connect(DATA/'fresh_AVI_exact_scores.sqlite')
    db.execute('CREATE TABLE IF NOT EXISTS scores (variant TEXT PRIMARY KEY, raw REAL, phred REAL, available INTEGER, source TEXT)')
    for key in wanted:
        row = old.execute('SELECT variant,raw,phred,available,source FROM scores WHERE variant=?', (key,)).fetchone()
        if row is not None:
            db.execute('INSERT OR IGNORE INTO scores VALUES (?,?,?,?,?)', row)
    old.close(); db.commit()
    present = {r[0] for r in db.execute('SELECT variant FROM scores')}
    groups = collections.defaultdict(set)
    for key in wanted-present:
        c, p, a, b = key.split(':'); p = int(p)
        groups[(c, (p-1)//16384)].add((p, a, b))
    local = threading.local()
    def query(item):
        (c, w), keys = item
        if not hasattr(local, 'tab'):
            local.tab = pysam.TabixFile(URL, index=str(index))
        found = {}
        for line in local.tab.fetch(c, min(k[0] for k in keys)-1, max(k[0] for k in keys)):
            f = line.split('\t'); k = (int(f[1]), f[2], f[3])
            if k in keys:
                value = (float(f[4]), float(f[5]))
                assert np.isfinite(value).all() and value[1] >= 0
                assert k not in found or found[k] == value
                found[k] = value
        return [(f'{c}:{p}:{a}:{b}', *found.get((p, a, b), (None, None)), int((p, a, b) in found), 'exact_public_AVI_V8_fresh_input_query') for p, a, b in keys]
    failures = []; done = 0
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        pending = {pool.submit(query, item): item[0] for item in groups.items()}
        for future in concurrent.futures.as_completed(pending):
            try:
                db.executemany('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)', future.result())
            except Exception as error:
                failures.append({'window': list(pending[future]), 'error': repr(error)})
            done += 1
            if done % 20 == 0:
                db.commit(); print('V8 FRESH AVI INPUTS', done, '/', len(groups), flush=True)
    db.commit()
    write_json(OUT/'fresh_AVI_query_attempt_receipt.json', {'windows': len(groups), 'failures': failures})
    assert not failures, failures
    rows = []
    for key in sorted(wanted):
        v = db.execute('SELECT raw,phred,available FROM scores WHERE variant=?', (key,)).fetchone()
        assert v is not None
        rows.append({'variant_GRCh38': key, 'raw': v[0] if v[2] else '', 'phred': v[1] if v[2] else '', 'available': v[2]})
    write_csv(DATA/'fresh_AVI_comparison_inputs.csv', rows)
    db.close()
    for name, digest in bound['outputs'].items():
        assert sha(V7/name) == digest
    write_json(OUT/'AVI_comparison_input_receipt.json', {
        'source': URL, 'index_MD5': INDEX_MD5,
        'designation': 'Permissive Use Downloadable AVI SNV composite only',
        'not_used': 'Detailed AlphaGenome API predictions, feature breakdowns or predicted labels as experimental truth',
        'requested_fresh_alleles': len(wanted), 'available': sum(r['available'] for r in rows),
        'new_query_windows': len(groups), 'maximum_parallel_queries': 2,
        'training_bank_source_sha256': sha(V7/'avi_bank.npy'),
        'V8_training_bank_sha256': sha(DATA/'avi_training_bank.npy'),
        'fresh_input_file_sha256': sha(DATA/'fresh_AVI_comparison_inputs.csv'),
        'prior_AVI_arrays_unchanged': True,
        'fresh_model_performance_calculated': False,
        'missing_values': 'Unavailable values retain missing flags; no zero biological-effect substitution.',
        'scope': 'Prepared inputs for a later matched with/without experiment. These values are a composite including other predictors; any gain cannot be assigned uniquely to AlphaGenome. The unchanged V7 decision still excludes AVI.'
    })
    print('V8 MATCHED AVI INPUTS READY', len(rows), flush=True)


if __name__ == '__main__':
    main()
