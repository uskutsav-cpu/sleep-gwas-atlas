#!/usr/bin/env python3
"""Repeatable read-only verification of sealed v4 arithmetic and review artifacts.

No estimator, raw GWAS read or project imports. SHA checks read processed gzip
bytes without decompression. Independent erfc bisection checks liability algebra.
Use --receipt to save a first verification receipt; repeated default runs print.
"""
import csv
import datetime
import hashlib
import json
import math
import sys
from decimal import Decimal
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(dict(check=name,pass_=True))


checks=[]
receipt=HERE/'qc_arithmetic_receipt_v4.json'
d=json.loads(receipt.read_text())
review=HERE.parent/'reviews/statistical_geneticist_critical_review_v4.json'
rd=json.loads(review.read_text())
for r in d['inputs']:
    check('arithmetic_input_identity_'+r['path'],sha(r['path'])==r['sha256_before']==r['sha256_after'])
for r in d['outputs']+rd['inputs']+rd['outputs']:
    check('artifact_identity_'+r['path'],sha(r['path'])==r['sha256'])
check('upstream_checks',len(d['checks'])==392 and all(r['pass_'] for r in d['checks']))
with (HERE/'liability_input_sensitivity_v4.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
check('two_real_target_traits',{r['trait'] for r in rows}=={'ms','melanoma'})
for r in rows:
    p,k=float(r['logged_case_fraction']),float(r['frozen_population_K'])
    # Independent standard-normal survival inversion, avoiding NormalDist.
    lo,hi=-12.,12.
    for _ in range(100):
        t=(lo+hi)/2
        if .5*math.erfc(t/math.sqrt(2))>k:
            lo=t
        else:
            hi=t
    t=(lo+hi)/2
    log_c=2*math.log(k)+2*math.log1p(-k)-math.log(p)-math.log1p(-p)+t*t+math.log(2*math.pi)
    c=math.exp(log_c)
    check('independent_erfc_conversion_'+r['trait'],math.isclose(c,float(r['original_conversion_factor']),rel_tol=4e-14))
    p_exact=Decimal(r['logged_case_fraction'])
    factor=Decimal(4)*p_exact*(Decimal(1)-p_exact)
    h=Decimal(r['original_reported_liability_h2'])*factor
    se=Decimal(r['original_reported_liability_SE'])*factor
    check('independent_decimal_h2_rescale_'+r['trait'],math.isclose(float(h),float(r['balanced_same_fit_liability_h2']),rel_tol=3e-15))
    check('independent_decimal_SE_rescale_'+r['trait'],math.isclose(float(se),float(r['balanced_same_fit_liability_SE']),rel_tol=3e-15))
    check('Z_invariance_'+r['trait'],math.isclose(float(h/se),float(r['unchanged_Z']),rel_tol=3e-15))
    check('actual_total_N_'+r['trait'],int(r['cases'])+int(r['controls'])==int(r['actual_total_N']))
result=dict(schema='v4_artifact_identity_and_independent_liability_arithmetic_verification',
            completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='PASS',
            verifier_sha256=sha(Path(__file__)),arithmetic_receipt_sha256=sha(receipt),review_receipt_sha256=sha(review),
            checks=checks,estimator_calls=0,dense_raw_reads=0,
            limitation='Verifies arithmetic and bound artifact identity, not native reproduction or scientific QC resolution.')
if sys.argv[1:]==['--receipt']:
    with (HERE/'qc_artifact_verification_v4.json').open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
elif sys.argv[1:]:
    raise SystemExit('Usage: verify_qc_artifacts_v4.py [--receipt]')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}|dict(passing_checks=len(checks)),indent=2))
