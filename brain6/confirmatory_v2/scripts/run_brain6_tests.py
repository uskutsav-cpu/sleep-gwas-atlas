#!/usr/bin/env python3
"""Run explicit source-free Brain6 subset; list all unexecuted integrations."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];A=ROOT/'brain6/confirmatory_v2'
m=json.loads((A/'qc/brain6_test_dependency_manifest.json').read_text())
cmd=[sys.executable,'-m','pytest','-q','brain6/scripts/tests','extensions/brain6/tests','brain6/confirmatory_v2/tests','-m','not native']
for t in m['integration_tests']:cmd+=['--deselect',t['collected_node_id']]
env=dict(os.environ,PYTHONPATH=str(ROOT/'extensions/brain6'))
r=subprocess.run(cmd,cwd=ROOT,env=env)
(A/'qc/brain6_test_partition.json').write_text(json.dumps({'command':cmd,'returncode':r.returncode,'data_dependent_tests_not_run':m['integration_tests'],'native_tests_not_run_in_source_free_suite':True,'scientific_interpretation':'Software validation only; original full integration suite remains required'},indent=2)+'\n')
sys.exit(r.returncode)
