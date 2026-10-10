"""Own tiny source-header identity edge: never takes production inputs."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import gzip

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
sys.path.insert(0,str(S))
ORIGINAL=P.parent/'scripts/01_harmonize.py';CANDIDATE=S/'core_bounded_harmonizer_v1.py'
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
before={str(x):sha(x) for x in [ORIGINAL,CANDIDATE,S/'core_column_reader_v3.py']}
assert before[str(CANDIDATE)]=='9aa49e2239c52d16cd206b15058000f7e0b5842d6103f920abd1013ebc1703d4'
root=R/'core_bounded_harmonizer_source_edge_controls_v1';root.mkdir(exist_ok=False)
config=root/'analysis_panel.tsv';config.write_text('trait_id\ttype\tbuild\tn_total\tsource_note\nx\tcontinuous\thg19\t1000\tOWN_SYNTHETIC_FIXTURE\n')
source=root/'source.tsv'
sealed='SNP\tCHR\tBP\tA1\tA2\tbeta\tb\tSE\tP\tN\nrs1\t1\t42\ta\tc\t1\t9\t1\t0.5\t1000\n'
unsealed='SNP\tCHR\tBP\tA1\tA2\tb\tSE\tP\tN\nrs1\t1\t42\ta\tc\t9\t1\t0.5\t1000\n'
source.write_text(sealed);expected=sha(source)
base=['--trait','x','--config',str(config),'--infile',str(source)]
spec=importlib.util.spec_from_file_location('bounded_candidate',CANDIDATE);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
original,_=c.load_original(ORIGINAL)
oldargv=sys.argv;sys.argv=[str(ORIGINAL),*base,'--outdir',str(root/'original')]
try:original.main()
finally:sys.argv=oldargv
source.write_text(unsealed)
loader=c.load_original;read=original.pd.read_csv;observed=[]
def load_patched(path):
    module,body=loader(path)
    def change_after_header(path,*args,**kwargs):
        result=read(path,*args,**kwargs)
        if str(path)==str(source) and kwargs.get('nrows')==0:
            observed.append(list(result.columns));source.write_text(sealed)
        return result
    module.pd.read_csv=change_after_header
    return module,body
c.load_original=load_patched
try:receipt=c.harmonize(ORIGINAL,[*base,'--outdir',str(root/'candidate')],root/'spool',1,expected)
finally:c.load_original=loader;original.pd.read_csv=read
oldbytes=gzip.decompress((root/'original/x.harmonized.tsv.gz').read_bytes())
newbytes=gzip.decompress((root/'candidate/x.harmonized.tsv.gz').read_bytes())
columns=json.loads((root/'spool/columns/column_preparation_manifest.json').read_text())
assert columns['source_sha256_before']==columns['source_sha256_after']==columns['expected_source_sha256']==expected==sha(source)
assert oldbytes!=newbytes and b'\t1\t1\t0.5\t1000\n' in oldbytes and b'\t9\t1\t0.5\t1000\n' in newbytes
assert observed and 'beta' not in observed[0]
after={p:sha(p) for p in before};assert before==after
result={'status':'MEASURED_FALSE_PASS_FULL_ADAPTER_PREFIX_HEADER_NOT_BOUND_TO_SEALED_RAW_IDENTITY','source_hash_expected':expected,'source_before_after_expected_all_match':True,'prefix_header_read_before_sealed_hash':observed,'candidate_selected_columns':columns['columns'],'original_decompressed_output':oldbytes.decode(),'candidate_decompressed_output':newbytes.decode(),'candidate_status':receipt['status'],'scientific_effect':'Sealed source has beta1 and redundant b9; original alias preference chooses beta. Adapter reads earlier unsealed header lacking beta then restores sealed bytes before freeze, chooses b9 and produces a different BETA with all later source hashes equal to expected.','minimal_fix':'Require expected source digest before any sniff/header/schema prefix read, retain final equality; owned immutable-source lifetime and metadata/dependency gates remain caller responsibilities.','source_bindings_before':before,'source_bindings_after':after,'production_body_reads':0,'production_mutex_or_workers_or_fits':0}
with (R/'core_bounded_harmonizer_source_edge_receipt_v1.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':result['status'],'prefix_header':observed,'candidate_columns':columns['columns']}))
