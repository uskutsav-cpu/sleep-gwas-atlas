#!/usr/bin/env python3
"""Integration check: a real LAVA input failure produces a row and nonzero exit."""
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V2 = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2')
RSCRIPT = Path('/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript')
RUNNER = ROOT / 'brain6/scripts/run_insomnia_native_n_pilot_v4.R'
CORE = ROOT / 'brain6/scripts/insomnia_native_n_pilot_v4_core.R'


def main() -> None:
    base = json.loads((V2 / 'worker_2.config.json').read_text())
    with tempfile.TemporaryDirectory(prefix='brain6-insomnia-v4-failure-') as tmp:
        root = Path(tmp)
        input_root = root / 'inputs'
        locus = input_root / 'locus_950'
        locus.mkdir(parents=True)
        (locus / 'input_info.tsv').write_text('phenotype\tcases\tcontrols\tfilename\ninsomnia\t109402\t277131\tinsomnia.sumstats.tsv.gz\n')
        # A correctly shaped nonempty shard with only one SNP must not be
        # special-cased as an empty locus; LAVA's own input gate should fail.
        with gzip.open(locus / 'insomnia.sumstats.tsv.gz', 'wt') as stream:
            stream.write('SNP\tA1\tA2\tZ\tN\nrs999999999\tA\tC\t1.2\t386533\n')
        base.update({'analysis_id': 'brain6_insomnia_native_n_trait_only_pilot_v4',
                     'scope': 'INTEGRATION_TEST_ONLY', 'locus_ids': ['950'],
                     'expected_loci': 1, 'input_rows': {'950': 1},
                     'input_root': str(input_root), 'core_path': str(CORE),
                     'output_tsv': str(root / 'output.tsv'),
                     'summary_json': str(root / 'summary.json')})
        config_path = root / 'config.json'
        config_path.write_text(json.dumps(base))
        env = dict(os.environ)
        env.update({'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'})
        result = subprocess.run([str(RSCRIPT), str(RUNNER), str(config_path)], cwd=ROOT,
                                env=env, capture_output=True, text=True, check=False)
        if result.returncode == 0:
            raise AssertionError(f'Genuine failure was swallowed: {result.stdout}')
        with (root / 'output.tsv').open(newline='') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        summary = json.loads((root / 'summary.json').read_text())
        if len(rows) != 1 or rows[0]['status'] != 'FAILED' or summary['failed'] != 1:
            raise AssertionError(f'Genuine failure not structured: {rows} / {summary}')
        print('INSOMNIA_V4_GENUINE_FAILURE_EXIT_TEST_PASS')


if __name__ == '__main__':
    main()
