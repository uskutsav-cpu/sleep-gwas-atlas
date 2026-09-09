#!/usr/bin/env python3
"""Measured synthetic I/O benchmark; NOT a real GWAS/native-method benchmark."""
from pathlib import Path
import argparse,json,math,platform,resource,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from brain6.io import require,write_json,write_tsv,sha256
from brain6.gwas import FIELDS,normalize,join_pair

p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--rows',type=int,default=100000)
a=p.parse_args();out=Path(a.out).resolve();require(not out.exists(),'Use a new output directory');out.mkdir(parents=True)
require(a.rows>=22,'At least 22 rows');start=time.perf_counter()
for k in [1,2]:
    def rows():
        for i in range(a.rows):
            z=((i%31)-15)/8;beta=z*.025;frequency=.3;a1,a2='A','C'
            if k==2 and i%2:a1,a2=a2,a1;beta=-beta;frequency=1-frequency
            yield dict(zip(FIELDS,[f'rs{i+1}',i%22+1,(i//22+1)*100,a1,a2,beta,.025,
                 math.erfc(abs(z)/math.sqrt(2)),100000,frequency,.99]))
    raw=out/f'raw{k}.tsv';write_tsv(raw,FIELDS,rows())
    card={'trait_id':f'synthetic{k}','study_id':'SYNTHETIC_BENCHMARK','phenotype_definition':'Artificial throughput fixture',
          'path':str(raw),'sha256':sha256(raw),'synthetic':True,'ancestry':'EUR','genome_build':'GRCh37',
          'n_semantics':'total','effect_scale':'beta','column_map':{x:x for x in FIELDS},
          'qc':{'require_eaf':True,'min_maf':.01,'min_n':100}}
    path=out/f'card{k}.json';write_json(path,card)
    normalize(path,out,f'normalized{k}',synthetic=True)
normalized=time.perf_counter()
join_pair(out/'normalized1',out/'normalized2',out,'pair',synthetic=True,min_variants=a.rows,min_overlap=1)
end=time.perf_counter();rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
result={'synthetic':True,'rows_per_trait':a.rows,'traits':2,'joined_rows':a.rows,
        'source_creation_and_normalization_seconds':normalized-start,'pair_join_seconds':end-normalized,
        'total_seconds':end-start,'maximum_resident_set_bytes':rss*(1 if sys.platform=='darwin' else 1024),
        'python':platform.python_version(),'platform':platform.platform(),
        'native_genetics_methods_executed':[],
        'scope':'Artificial parsing, normalization, SQLite indexing/join, allele alignment and hashing only. Not a whole-genome or R-method performance prediction.'}
write_json(out/'benchmark.json',result);print(json.dumps(result,indent=2))
