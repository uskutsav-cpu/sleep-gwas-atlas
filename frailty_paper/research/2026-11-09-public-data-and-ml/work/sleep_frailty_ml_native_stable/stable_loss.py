"""Same measured-data objective with fixed-shape masked reductions."""
import torch
import torch.nn.functional as F


def masked_mean(value,mask):
    return (value*mask.to(value.dtype)).sum()/mask.sum().clamp_min(1).to(value.dtype)


def activity_loss(prediction,batch,weights,quality,nmain,nlocal,ncross):
    p={k:v[:nmain] for k,v in prediction.items()};w=weights[:nmain];q=quality[:nmain]
    result=(w*abs(p['magnitude']-batch['effect'][:nmain])).mean()
    signed=torch.isfinite(batch['signed_effect'][:nmain])
    target=torch.where(signed,batch['signed_effect'][:nmain],torch.zeros_like(batch['signed_effect'][:nmain]))
    result=result+.25*masked_mean(w*q*abs(p['signed_effect']-target),signed)
    active=torch.isfinite(batch['activity_reference'][:nmain])&torch.isfinite(batch['activity_alternate'][:nmain])
    a=torch.where(active,batch['activity_reference'][:nmain]-batch['activity_center'][:nmain],torch.zeros_like(batch['activity_reference'][:nmain]))
    b=torch.where(active,batch['activity_alternate'][:nmain]-batch['activity_center'][:nmain],torch.zeros_like(batch['activity_alternate'][:nmain]))
    ea=F.smooth_l1_loss(p['activity_reference'],a,reduction='none');eb=F.smooth_l1_loss(p['activity_alternate'],b,reduction='none')
    result=result+.05*masked_mean(w*(ea+eb)*.5,active)
    result=result+.05*masked_mean(w*q*F.smooth_l1_loss(p['signed_effect'],b-a,reduction='none'),active)
    for count,offset in [(nlocal,nmain),(ncross,nmain+2*nlocal)]:
        if count:
            pp=prediction['magnitude'][offset:offset+2*count].reshape(count,2);y=batch['effect'][offset:offset+2*count].reshape(count,2)
            qq=quality[offset:offset+2*count].reshape(count,2).prod(1).sqrt()
            result=result+.2*(qq*abs((pp[:,0]-pp[:,1])-(y[:,0]-y[:,1]))).mean()
    return result
