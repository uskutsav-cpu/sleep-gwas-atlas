#!/usr/bin/env python3
"""Result-free stable-dict rsID equivalence and bounded Python3.9 timings."""
import contextlib
import hashlib
import importlib.util
import io
import json
import random
from pathlib import Path
import resource
import statistics
import sys
import time

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
sys.path.insert(0,str(S))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,p):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def save(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')

def main():
    if not sys.version.startswith('3.9.23'):raise RuntimeError('DECLARED_PYTHON3_9_REQUIRED')
    started=time.monotonic();oldpath=S/'56_preprocess_finngen_insomnia_feasibility_v2.py';newpath=S/'56_preprocess_finngen_insomnia_feasibility_v3.py'
    old,new=module('_old_rsid',oldpath),module('_new_rsid',newpath)
    code_before={str(x):sha(x) for x in [oldpath,newpath,S/'57_verify_finngen_row_controls_v3.py']}
    equivalent_source=newpath.read_text()==oldpath.read_text().replace('set(RSID.findall(value)).intersection(hm3)','{rsid for rsid in RSID.findall(value) if rsid in hm3}')
    lookup={('rs'+str(i)):None for i in range(30)}
    values=['','NA','rs1','rs999999','rs1;rs1','rs1,rs2','rs1x','xrs1','RS1','rs1/rs999999','rs0001','1rs1','rs1_foo','rs1.rs2']
    rng=random.Random(7307)
    for _ in range(500):
        tokens=[rng.choice(['rs'+str(i) for i in range(50)]+['NA','rs1x','xrs1','RS1','']) for j in range(rng.randrange(9))]
        values.append(rng.choice([';',',',' ','|','/']).join(tokens))
    equivalence=[]
    for value in values:
        a=old.unique_hm3_rsid(value,lookup);b=new.unique_hm3_rsid(value,lookup)
        equivalence.append(dict(value=value,old=a,new=b,passed=a==b))
    timings=[]
    for size in [1000,100000,400000]:
        hm3={'rs'+str(i):None for i in range(size)}
        for value in ['', 'NA', 'rs9999999','rs1','rs1;rs1','rs1,rs2']:
            samples={}
            for name,fn in [('old',old.unique_hm3_rsid),('new',new.unique_hm3_rsid)]:
                t=[]
                for _ in range(3):
                    begin=time.perf_counter()
                    for j in range(20):fn(value,hm3)
                    t.append((time.perf_counter()-begin)/20)
                samples[name]=statistics.median(t)
            timings.append(dict(dictionary_keys=size,value=value,median_seconds_per_call=samples,
                                old_over_new_ratio=samples['old']/samples['new'],results_equal=old.unique_hm3_rsid(value,hm3)==new.unique_hm3_rsid(value,hm3)))
        del hm3
    controls=module('_original34_v3',S/'57_verify_finngen_row_controls_v3.py')
    row_receipt=R/'independent_finngen_original34_row_recheck_v4_3_7.json'
    def redirected(path,record):
        if Path(path)!=P/'logs/finngen_row_controls_v4_3_7.json':raise RuntimeError('UNEXPECTED_ORIGINAL_CONTROL_OUTPUT')
        save(row_receipt,record)
    controls.write_new=redirected
    with contextlib.redirect_stdout(io.StringIO()):controls.main()
    rerun=json.loads(row_receipt.read_text());frozen=json.loads((P/'logs/finngen_row_controls_v4_3_7.json').read_text())
    same34=rerun['checks']==frozen['checks'] and len(rerun['checks'])==34 and all(x['passed'] for x in rerun['checks'])
    result=dict(schema='independent_stable_dict_rsid_equivalence_python39_v4_3_7',python=sys.version,
        python_executable=sys.executable,exact_only_set_expression_source_change=equivalent_source,
        equality_checks=equivalence,all514_equality_checks_pass=all(x['passed'] for x in equivalence),
        timings=timings,all_timing_cases_equal=all(x['results_equal'] for x in timings),
        original34_row_controls_rerun_identical=same34,row_recheck_sha256=sha(row_receipt),
        frozen_original_row_receipt_sha256=sha(P/'logs/finngen_row_controls_v4_3_7.json'),
        reviewed_code_sha256=code_before,reviewed_code_unchanged=all(sha(k)==v for k,v in code_before.items()),
        max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,elapsed_seconds=time.monotonic()-started,
        semantics_qualification='Stable ordinary dict with string keys, no concurrent mutation; dictionary values irrelevant. Not a real source-content or completed-stream replay.',
        real_body_reads=False,actual_scientific_worker_or_fit=False,network=False,heavy_mutex=False)
    save(R/'independent_finngen_rsid_equivalence_receipt_v4_3_7.json',result)
    print(json.dumps(dict(equivalence=all(x['passed'] for x in equivalence),row34_identical=same34,
        peak_RSS_bytes=result['max_RSS_bytes'],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
