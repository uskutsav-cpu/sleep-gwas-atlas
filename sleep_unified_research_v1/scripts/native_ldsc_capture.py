#!/usr/bin/env python3
"""Run pinned LDSC and capture full precision without changing its estimator.

Instrumentation wraps return values only. Stock logs and block-delete values
are retained. Input paths and their hash ledger must be admitted by the runner.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import runpy
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--ldsc-dir', type=Path, required=True)
known, forwarded = parser.parse_known_args()
code = known.ldsc_dir.resolve()
sys.path.insert(0, str(code))
import ldscore.sumstats as sumstats
import numpy as np
import pandas as pd
import scipy

def clean(value):
    if isinstance(value, (np.integer, np.floating)): value = value.item()
    if isinstance(value, float) and not np.isfinite(value): return None
    if isinstance(value, np.ndarray): return value.tolist()
    return value

def attrs(obj, fields):
    return {f: clean(getattr(obj, f, None)) for f in fields}

original_rg, original_h2 = sumstats.estimate_rg, sumstats.estimate_h2

def capture_rg(args, log):
    results = original_rg(args, log)
    paths = args.rg.split(',')
    rows = []
    for path, obj in zip(paths[1:], results):
        row = {'p1': paths[0], 'p2': path}
        if obj is None:
            row['status'] = 'ESTIMATION_FAILED'
        else:
            row.update(attrs(obj, ['rg_ratio','rg_jknife','rg_se','z','p']))
            for prefix, value in [('hsq1',obj.hsq1),('hsq2',obj.hsq2),('gencov',obj.gencov)]:
                row[prefix] = attrs(value,['tot','tot_se','intercept','intercept_se','mean_chisq','lambda_gc','ratio','ratio_se'])
            row['status'] = 'NATIVE_ESTIMATE_RETURNED'
        rows.append(row)
    write_receipt(args, rows)
    return results

def capture_h2(args, log):
    result = original_h2(args, log)
    # estimate_h2 returns the fitted regression object at the pinned revision.
    write_receipt(args, [dict(input=args.h2, **attrs(result,['tot','tot_se','intercept','intercept_se','mean_chisq','lambda_gc','ratio','ratio_se']))])
    return result

def write_receipt(args, rows):
    result = {'completed_utc':datetime.now(timezone.utc).isoformat(),
              'instrumentation':'return-value capture; numerical estimator unmodified',
              'ldsc_dir':str(code), 'python':sys.version,
              'libraries':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
              'arguments':vars(args), 'estimates':rows}
    Path(args.out+'.full_precision.json').write_text(json.dumps(result,indent=2,default=str,allow_nan=False)+'\n')

sumstats.estimate_rg = capture_rg
sumstats.estimate_h2 = capture_h2
sys.argv = [str(code/'ldsc.py')] + forwarded
runpy.run_path(str(code/'ldsc.py'),run_name='__main__')
