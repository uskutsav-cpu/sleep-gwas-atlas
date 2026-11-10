"""Verify the distributed snapshot and optionally replay its saved predictions."""
import argparse
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
NATIVE = Path('outputs/sleep_frailty_ml_resumable')


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024**2), b''):
            value.update(block)
    return value.hexdigest()


def verify():
    rows = list(csv.DictReader((ROOT / 'release_manifest.tsv').open(), delimiter='\t'))
    assert rows and len({row['path'] for row in rows}) == len(rows)
    for row in rows:
        relative = Path(row['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path = ROOT / relative
        assert path.is_file() and not path.is_symlink(), row['path']
        assert path.stat().st_size == int(row['bytes']), row['path']
        assert digest(path) == row['sha256'], row['path']
    core = ROOT / 'work/sleep_frailty_ml_resumable'
    scripts = list(core.rglob('*.py'))
    for script in scripts:
        ast.parse(script.read_text(), filename=str(script))
    decision = json.loads((ROOT / NATIVE / 'native_comparison_independent_decision_audit.json').read_text())
    assert decision['passed'] and not decision['scientific_gate_passed']
    assert decision['all_eight_fits_complete_and_matched']
    assert decision['all_24_saved_magnitude_results_independently_checked']
    confidence = json.loads((ROOT / 'outputs/sleep_frailty_ml_confidence_revision/confidence_method_decision.json').read_text())
    assert confidence['qualifying_sources'] == 0 and not confidence['new_models_promoted']
    print(json.dumps({'snapshot_hashes_verified': len(rows), 'current_core_AST_valid': len(scripts),
                      'native_accuracy_gate_failed': True, 'confidence_candidate_promoted': False}), flush=True)


def replay():
    # Auditors write receipts; isolate those writes from the distributed snapshot.
    with tempfile.TemporaryDirectory(prefix='sleep-frailty-replay-') as directory:
        target = Path(directory)
        source_scripts = Path('work/sleep_frailty_ml_resumable')
        for name in ['audit_saved_effect_size_metrics.py', 'audit_native_comparison_decision.py']:
            destination = target / source_scripts / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / source_scripts / name, destination)
        data = Path('outputs/sleep_frailty_ml_next/coverage_training_arrays.npz')
        (target / data).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / data, target / data)
        (target / NATIVE / 'development').mkdir(parents=True)
        for path in (ROOT / NATIVE).glob('*.json'):
            shutil.copyfile(path, target / NATIVE / path.name)
        for path in (ROOT / NATIVE / 'development').iterdir():
            if path.suffix in {'.json', '.npz', '.safetensors'}:
                shutil.copyfile(path, target / NATIVE / 'development' / path.name)
        env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'OMP_NUM_THREADS': '2',
               'OPENBLAS_NUM_THREADS': '2', 'VECLIB_MAXIMUM_THREADS': '2'}
        count = 0
        for source in ['none', 'Salomon', 'Kosicki2025', 'Retallick2026']:
            for recipe in ['frozen_reader', 'adapt_last_stage']:
                for epoch in (1, 2, 3):
                    result = subprocess.run(
                        [sys.executable, str(target / source_scripts / 'audit_saved_effect_size_metrics.py'),
                         '--recipe', recipe, '--source', source, '--epoch', str(epoch)],
                        env=env, check=True, capture_output=True, text=True)
                    assert json.loads(result.stdout)['passed']
                    count += 1
        subprocess.run([sys.executable, str(target / source_scripts / 'audit_native_comparison_decision.py')],
                       env=env, check=True, stdout=subprocess.DEVNULL)
        print(json.dumps({'actual_saved_magnitude_replays_passed': count,
                          'locked_stopping_decision_replayed': True,
                          'reserved_outcomes_used': False,
                          'model_training_performed': False}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay-native', action='store_true')
    arguments = parser.parse_args()
    verify()
    if arguments.replay_native:
        replay()
