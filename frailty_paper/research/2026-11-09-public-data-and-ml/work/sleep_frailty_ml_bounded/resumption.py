"""Atomic, weights-only-safe recovery of training state and randomness."""
import json,os
from pathlib import Path
import torch


def cpu_tree(value):
    if isinstance(value,torch.Tensor):return value.detach().cpu().clone()
    if isinstance(value,dict):return {k:cpu_tree(v) for k,v in value.items()}
    if isinstance(value,list):return [cpu_tree(v) for v in value]
    if isinstance(value,tuple):return tuple(cpu_tree(v) for v in value)
    return value


def save_recovery(path,model,optimizer,device,identity,progress):
    path=Path(path);path.parent.mkdir(exist_ok=True)
    if device=='mps':torch.mps.synchronize()
    parameters={k:p.detach().cpu().clone() for k,p in model.named_parameters() if p.requires_grad}
    payload={'format_version':1,'identity_json':json.dumps(identity,sort_keys=True),
        'progress_json':json.dumps(progress,sort_keys=True),'parameters':parameters,
        'optimizer':cpu_tree(optimizer.state_dict()),'CPU_rng':torch.get_rng_state(),
        'MPS_rng':torch.mps.get_rng_state() if device=='mps' else torch.empty(0,dtype=torch.uint8)}
    temporary=path.with_name(path.name+'.pending')
    with temporary.open('wb') as handle:
        torch.save(payload,handle);handle.flush();os.fsync(handle.fileno())
    os.replace(temporary,path)


def load_recovery(path,model,optimizer,device,identity):
    # Only this project's own primitive/tensor payload is accepted. No arbitrary
    # Python objects or unrestricted pickle deserialization is used.
    payload=torch.load(path,map_location='cpu',weights_only=True)
    assert payload['format_version']==1
    assert json.loads(payload['identity_json'])==identity,'Recovery identity/source mismatch'
    expected={k for k,p in model.named_parameters() if p.requires_grad}
    assert set(payload['parameters'])==expected
    actual=dict(model.named_parameters())
    with torch.no_grad():
        for key,value in payload['parameters'].items():actual[key].copy_(value.to(actual[key].device))
    optimizer.load_state_dict(payload['optimizer'])
    torch.set_rng_state(payload['CPU_rng'])
    if device=='mps':torch.mps.set_rng_state(payload['MPS_rng'])
    progress=json.loads(payload['progress_json'])
    assert progress['updates']>=0
    return progress
