"""Pinned author architecture and tensor-only extraction of its real weights."""
import importlib.util, pickletools, zipfile
from v8_common import *
import torch
from safetensors.torch import load_file, save_file

torch.set_num_threads(2)


def extract_weights():
    source = SOURCES/'legnet_models_example.cpkt'
    target = MODELS/'author_WTC11_original.safetensors'
    receipt_path = OUT/'pretrained_weight_extraction_audit.json'
    if target.exists():
        receipt = json.loads(receipt_path.read_text())
        assert sha(source) == receipt['original_checkpoint_sha256']
        assert sha(target) == receipt['tensor_file_sha256']
        return target
    # The author's training checkpoint also serializes optimizer scheduler
    # objects. Do not execute them or disable PyTorch's weights-only loader.
    # Its state dictionary is complete immediately before the top-level
    # "loops" field. STOP at that exact boundary returns the completed tensor
    # dictionary, retaining all original tensor storage bytes unchanged.
    tensor_archive = MODELS/'author_state_only_extraction.pt'
    with zipfile.ZipFile(source) as z:
        assert z.testzip() is None
        name = next(n for n in z.namelist() if n.endswith('/data.pkl'))
        data = z.read(name)
        ops = list(pickletools.genops(data))
        boundaries = [i for i, (op, arg, pos) in enumerate(ops) if op.name == 'BINUNICODE' and arg == 'loops']
        assert len(boundaries) == 1
        i = boundaries[0]
        assert ops[i-1][0].name == 'BUILD'
        boundary = ops[i][2]
        trimmed = data[:boundary]+b'.'
        globals_before = [arg for op, arg, pos in pickletools.genops(trimmed) if op.name == 'GLOBAL']
        assert set(globals_before) <= {'collections OrderedDict', 'torch._utils _rebuild_tensor_v2', 'torch FloatStorage', 'torch LongStorage'}
        with zipfile.ZipFile(tensor_archive, 'w', compression=zipfile.ZIP_STORED) as dest:
            for member in z.infolist():
                dest.writestr(member.filename, trimmed if member.filename == name else z.read(member.filename))
    state = torch.load(tensor_archive, map_location='cpu', weights_only=True)
    assert state and all(isinstance(k, str) and k.startswith('model.') and isinstance(v, torch.Tensor) for k, v in state.items())
    with zipfile.ZipFile(source) as old, zipfile.ZipFile(tensor_archive) as new:
        storages = [n for n in old.namelist() if '/data/' in n]
        assert all(hashlib.sha256(old.read(n)).digest() == hashlib.sha256(new.read(n)).digest() for n in storages)
    state = {k.removeprefix('model.'): v.contiguous() for k, v in state.items()}
    save_file(state, target)
    recovered = load_file(target)
    assert set(recovered) == set(state)
    assert all(torch.equal(recovered[k], v) for k, v in state.items())
    write_json(receipt_path, {
        'original_checkpoint': str(source), 'original_checkpoint_sha256': sha(source),
        'tensor_file': str(target), 'tensor_file_sha256': sha(target),
        'tensor_count': len(state), 'all_serialized_tensor_storage_bytes_preserved': True,
        'all_reloaded_tensors_exactly_equal': True, 'optimizer_or_scheduler_objects_deserialized': False,
        'unsafe_pickle_loader_enabled': False, 'static_prefix_boundary': boundary,
        'allowed_globals_in_state_prefix': globals_before,
        'author_checkpoint_epoch_metadata': 22,
        'exact_pretraining_fold_unknown': True,
        'scope': 'Real author-distributed WTC11 pretrained tensor values, not newly initialized or third-party reconstructed weights. This numerical extraction does not establish predictive accuracy or exclusion of pretraining sequence overlap.'
    })
    return target


def base_model():
    extract_weights()
    path = SOURCES/'legnet_model.py'
    spec = importlib.util.spec_from_file_location('pinned_author_legnet', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cfg = json.loads((SOURCES/'legnet_models_example.cfg').read_text())
    keys = ['stem_ch', 'stem_ks', 'ef_ks', 'ef_block_sizes', 'pool_sizes', 'resize_factor']
    model = module.LegNet(in_ch=4, **{k: cfg[k] for k in keys})
    model.load_state_dict(load_file(MODELS/'author_WTC11_original.safetensors'), strict=True)
    return model.eval()


def encode(sequences, device='cpu', dtype=torch.float32):
    # The author uses A,G,C,T channel order, NOT the common A,C,G,T order.
    assert len({len(s) for s in sequences}) == 1
    lut = np.full(256, 255, dtype='uint8')
    lut[[ord(c) for c in 'AGCT']] = [0, 1, 2, 3]
    tokens = np.stack([lut[np.frombuffer(s.encode(), dtype='uint8')] for s in sequences])
    assert np.max(tokens) < 4
    return torch.nn.functional.one_hot(torch.from_numpy(tokens.astype('int64')), 4).permute(0, 2, 1).to(device=device, dtype=dtype)


def crop_pair(row, width=150):
    a, b = row['reference'], row['alternate']
    p = int(row['changed_position'])
    assert len(a) == len(b) and len(a) >= width
    start = max(0, min(p-width//2, len(a)-width))
    a, b = a[start:start+width], b[start:start+width]
    assert sum(x != y for x, y in zip(a, b)) == 1
    return a, b
