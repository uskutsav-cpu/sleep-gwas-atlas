#!/usr/bin/env python3
"""Separate explicit archived-data integration tests from source-free CI.

Never counts unavailable integration tests as passing. Original assertions and
test files are unchanged. Any unclassified failure remains a CI failure.
"""
import argparse,json,sys,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
AREA=ROOT/'brain6/confirmatory_v2'

def flatten(suite):
    for test in suite:
        if isinstance(test,unittest.TestSuite):yield from flatten(test)
        else:yield test

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--suite',choices=['source-free','integration'],default='source-free')
    parser.add_argument('--report',type=Path,default=AREA/'qc/contract_test_partition.json');args=parser.parse_args()
    specs=json.loads((AREA/'qc/test_dependency_manifest.json').read_text())['integration_tests']
    manifest={r['test_id']:r for r in specs}
    loader=unittest.TestLoader();all_tests=list(flatten(loader.discover(str(ROOT/'tests'))))
    integration=[];unit=[]
    for t in all_tests:
        ident=t.id();module=ident.split('.')[0]
        if ident.startswith('unittest.loader._FailedTest.'):module=ident.split('.')[-1]
        key=ident if ident in manifest else module+'.*'
        (integration if key in manifest else unit).append(t)
    report={'mode':args.suite,'source_free_test_count':len(unit),'integration_test_count':len(integration),
        'integration_status':'NOT_RUN_ARCHIVED_DATA_OR_RUNTIME_REQUIRED' if args.suite=='source-free' else 'RUN_REQUESTED',
        'integration_tests':[t.id() for t in integration],'assertions_removed':0,
        'interpretation':'A source-free pass does not validate production statistical results'}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    selected=unit if args.suite=='source-free' else integration
    result=unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(selected))
    report.update(tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
                  existing_skips=[{'test':t.id(),'reason':reason} for t,reason in result.skipped],success=result.wasSuccessful())
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'Explicit partition: {len(unit)} source-free; {len(integration)} archived-data/native integration tests. Integration is not a source-free success.')
    return 0 if result.wasSuccessful() else 1

if __name__=='__main__':sys.exit(main())
