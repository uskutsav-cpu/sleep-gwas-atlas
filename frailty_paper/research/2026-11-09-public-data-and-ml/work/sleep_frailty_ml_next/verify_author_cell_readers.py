"""Recover original tensor values without deserializing optimizer/scheduler objects."""
import io,json,hashlib,zipfile,pickletools,importlib.util
from pathlib import Path
import numpy as np
import torch
from safetensors.torch import load_file,save_file

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'
ARCH=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v8/sources/legnet_model.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    torch.set_num_threads(2)
    receipt=json.loads((OUT/'cell_matched_pretrained_retrieval_receipt.json').read_text())
    for f in receipt['files']:assert sha(f['path'])==f['sha256']
    assert sha(ARCH)=='8a5650b788662db07fda0e82566697b0aa1b81a1e6ff7a53448f138d00c40975'
    spec=importlib.util.spec_from_file_location('pinned_author_architecture_for_new_readers',ARCH);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    checks={}
    for cell in ['K562','HepG2','WTC11']:
        source=OUT/f'{cell}_author_test1_val2.ckpt';target=OUT/f'{cell}_author_test1_val2.safetensors';memory=io.BytesIO()
        with zipfile.ZipFile(source) as old:
            assert old.testzip() is None;name=next(n for n in old.namelist() if n.endswith('/data.pkl'));raw=old.read(name);ops=list(pickletools.genops(raw))
            locations=[i for i,(op,arg,pos) in enumerate(ops) if op.name=='BINUNICODE' and arg=='loops'];assert len(locations)==1;index=locations[0];assert ops[index-1][0].name=='BUILD';boundary=ops[index][2]
            prefix=raw[:boundary]+b'.';globals_before={arg for op,arg,_ in pickletools.genops(prefix) if op.name=='GLOBAL'}
            assert globals_before<={'collections OrderedDict','torch._utils _rebuild_tensor_v2','torch FloatStorage','torch LongStorage'}
            with zipfile.ZipFile(memory,'w',compression=zipfile.ZIP_STORED) as clean:
                for member in old.infolist():clean.writestr(member.filename,prefix if member.filename==name else old.read(member.filename))
            memory.seek(0)
            with zipfile.ZipFile(memory) as clean:
                storages=[n for n in old.namelist() if '/data/' in n]
                assert all(hashlib.sha256(old.read(n)).digest()==hashlib.sha256(clean.read(n)).digest() for n in storages)
        memory.seek(0);original=torch.load(memory,map_location='cpu',weights_only=True)
        assert original and all(k.startswith('model.') and isinstance(v,torch.Tensor) for k,v in original.items())
        state={k.removeprefix('model.'):v.contiguous() for k,v in original.items()}
        if not target.exists():save_file(state,target)
        reloaded=load_file(target);assert set(reloaded)==set(state) and all(torch.equal(reloaded[k],state[k]) for k in state)
        cfg=json.loads((OUT/f'{cell}_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor'];model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(state,strict=True);model.eval()
        # Synthetic unambiguous DNA is only a tensor/shape check. No new
        # experimental DNA or effects are used to choose a reader here.
        torch.manual_seed(20261108);dimensions={}
        with torch.no_grad():
            for length in [148,150,200,230]:
                tokens=torch.randint(4,(4,length));x=torch.nn.functional.one_hot(tokens,4).permute(0,2,1).float();y=model(x);assert y.shape==(4,) and torch.isfinite(y).all();dimensions[str(length)]='finite four-case forward pass'
        checks[cell]={'checkpoint_sha256':sha(source),'tensor_file_sha256':sha(target),'tensor_count':len(state),'architecture_parameter_count':sum(p.numel() for p in model.parameters()),'author_test_fold':1,'author_validation_fold':2,
                      'all_original_tensor_storage_bytes_unchanged':True,'all_loaded_tensors_equal_original_values':True,'unsafe_optimizer_or_scheduler_deserialization':False,'permitted_prefix_globals':sorted(globals_before),'prefix_boundary':boundary,'synthetic_input_checks':dimensions}
    (OUT/'cell_reader_tensor_and_architecture_audit.json').write_text(json.dumps({'passed':True,'author_architecture_sha256':sha(ARCH),'checks':checks,'source_ensemble_scope':'One fixed author fold per cell, not the full 90-model ensemble.','new_experimental_cases_scored':False,'stronger_prediction_established':False,'future_pretraining_exposure_audit_required':True},indent=2)+'\n')
    print('THREE REAL CELL-MATCHED AUTHOR READERS VERIFIED; PREDICTIVE VALUE STILL UNTESTED',flush=True)


if __name__=='__main__':main()
