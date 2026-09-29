#!/usr/bin/env python3
"""Read-only structural comparison of historical and current Track B input gates."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ARCHIVE = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
OLD_GATE = ARCHIVE / 'results/track_b/pleiotropy/input_gate.lock.json'
NEW_GATE = ROOT / 'results/track_b/pleiotropy/input_gate.lock.json'
OLD_READY = ARCHIVE / 'results/track_b/pleiotropy/readiness_gate.tsv'
NEW_READY = ROOT / 'results/track_b/pleiotropy/readiness_gate.tsv'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def tsv(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    return {x['component_id']:x for x in rows}

def write_once(path,data):
    if path.exists():
        assert path.read_text()==data
    else:
        path.write_text(data)

def main():
    old=json.loads(OLD_GATE.read_text()); new=json.loads(NEW_GATE.read_text())
    all_fields=sorted(set(old)|set(new))
    gate_delta={k:{'historical':old.get(k),'current':new.get(k)} for k in all_fields if old.get(k)!=new.get(k)}
    assert set(gate_delta)=={'readiness_gate_sha256'}
    assert old['readiness_gate_sha256']==sha(OLD_READY)
    assert new['readiness_gate_sha256']==sha(NEW_READY)
    a,b=tsv(OLD_READY),tsv(NEW_READY)
    assert set(a)==set(b)
    ready_delta=[]
    for component in sorted(a):
        assert set(a[component])==set(b[component])
        for field in a[component]:
            if a[component][field]!=b[component][field]:
                ready_delta.append({'component_id':component,'field':field,
                                    'historical':a[component][field], 'current':b[component][field]})
    assert len(ready_delta)==2
    assert {x['component_id'] for x in ready_delta}=={'PLACO_PLUS','METHOD_UNION_08_09'}
    assert all(x['field']=='blockers' and x['current']==x['historical']+';ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION' for x in ready_delta)
    result={'analysis_id':'brain6-track-b-gate-delta-v1','status':'READ_ONLY_DIAGNOSTIC_NOT_ADMISSION',
            'historical_input_gate_sha256':sha(OLD_GATE),'current_input_gate_sha256':sha(NEW_GATE),
            'historical_readiness_sha256':sha(OLD_READY),'current_readiness_sha256':sha(NEW_READY),
            'gate_delta':gate_delta,'readiness_delta':ready_delta,
            'interpretation':'The current input-gate identity changed only through two appended readiness blocker strings; this does not itself prove protected-v3 admission.'}
    out=OUT/'gate_delta_audit.json'
    write_once(out,json.dumps(result,indent=2,sort_keys=True)+'\n')
    report='''# Historical-to-current Track B gate delta

The two input_gate.lock.json objects have identical fields and values except readiness_gate_sha256. Both readiness hashes verify against their respective TSV bytes. The TSVs have the same component rows and columns; exactly two blockers cells differ. In both, current state appends ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION to the historical string. All other fields, including the five dense input identities and contract lock, are identical.

This narrows the v1 gate mismatch to readiness blocker text, but it does not admit historical B to frozen Brain6 v3. Any exception for this historical gate lineage still requires the proposed, explicitly authorized admission amendment and independent validator. See gate_delta_audit.json for exact before/after values and hashes.
'''
    write_once(OUT/'gate_delta_report.md',report)
    print(json.dumps({'input_gate_differences':len(gate_delta),'readiness_cell_differences':len(ready_delta)}))

if __name__=='__main__':main()
