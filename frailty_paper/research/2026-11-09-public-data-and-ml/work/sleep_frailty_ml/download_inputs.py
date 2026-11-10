"""Fetch only pinned, public source inputs and verify their recorded hashes."""
from pathlib import Path
import hashlib,subprocess
WORK=Path(__file__).resolve().parent;ROOT=WORK.parents[1]
sources=[
 ('https://raw.githubusercontent.com/kircherlab/80K-Analysis/a686f1d552f97ebcc0224057b5a835002f7a3dc8/variant_based_analysis/data/80k_tested_variant_annotation_table_with_readout.tsv.gz',ROOT/'work/salomon_neuronal_variant_readout.tsv.gz','fce979affffcc276c3e2c03bb4d4cb9307af7dbb5b6807be68af91a90bb15f48'),
 ('https://huggingface.co/LongSafari/hyenadna-tiny-1k-seqlen-hf/resolve/e8c1effa8673814e257e627d2e1eda9ea5a373f6/model.safetensors',WORK/'hyenadna/model.safetensors','5ce2146c21e9c4baa6bddc4998fd3d029903ae84a563bf80218644082194a12d')]
for url,path,sha in sources:
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:raise ValueError(f'Existing source checksum conflict: {path}')
        print('Verified existing',path);continue
    path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.download')
    subprocess.run(['curl','--fail','--location','--show-error','--max-time','180','--output',str(temp),url],check=True)
    if hashlib.sha256(temp.read_bytes()).hexdigest()!=sha:raise ValueError(f'Downloaded source checksum conflict: {temp}')
    temp.replace(path);print('Downloaded and verified',path)
