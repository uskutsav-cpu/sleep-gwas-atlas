"""Independent sequence arithmetic, transformations and magnitude-backprop checks."""
import gzip,collections,itertools
from common import *
from sequence_features import Encoder,tokenize,patterns,model,motif_model,mx
from network import initialize,forward,loss
def main():
    receipt=json.loads((OUT/'feature_receipt.json').read_text());assert sha(DATA/'feature_bank.npy')==receipt['sha256']
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    ids=[]
    for kind in ['v6','kosicki','hepg2','new']:ids.extend([i for i,r in enumerate(bank) if r['context_feature_source']==kind][:3])
    encoder=Encoder();x=np.load(DATA/'feature_bank.npy',mmap_mode='r');checks=[]
    for i in ids:
        r=bank[i];p=int(r['changed_position']);start=max(0,min(p-25,len(r['reference'])-51));a=tokenize(r['reference'][start:start+51])[None,:].astype('int32');b=tokenize(r['alternate'][start:start+51])[None,:].astype('int32');pp=np.array([p-start],dtype='int32');expected=np.asarray(model.strand_features(a,b,pp,encoder.numpy_weights))[:,128:];m=motif_model.scan_numpy(a,b,pp,encoder.numpy_motifs)
        width=min(151,len(r['reference']));start=max(0,min(p-75,len(r['reference'])-width));aa=tokenize(r['reference'][start:start+width])[None,:].astype('int32');bb=tokenize(r['alternate'][start:start+width])[None,:].astype('int32');context=np.asarray(model.strand_features(aa,bb,np.array([p-start]),encoder.numpy_weights))[:,:128]
        errors={'variant51_error':float(np.max(abs(expected-x[i:i+1,128:512]))),'motif51_error':float(np.max(abs(m-x[i:i+1,680:]))),'context_error':float(np.max(abs(context-x[i:i+1,:128])))}
        def literal(tokens):
            seq=''.join('ACGT'[z-7] for z in tokens);result=[]
            for k in [1,2,3]:
                count=collections.Counter(seq[j:j+k] for j in range(len(seq)-k+1));result.extend(count[''.join(w)] for w in itertools.product('ACGT',repeat=k))
            return np.array(result)
        den=np.r_[np.full(4,51),np.full(16,50),np.full(64,49)]
        def vector(u,v):
            au,bu=literal(u),literal(v);return np.r_[.5*(au+bu)/den,abs(au-bu)]
        pat=.5*(vector(a[0],b[0])+vector(17-a[0,::-1],17-b[0,::-1]));errors['literal_pattern_error']=float(np.max(abs(pat-x[i,512:680])));assert max(errors.values())<5e-5,(i,errors);checks.append({'bank_id':i,'source_context':r['context_feature_source'],**errors})
    # Same DNA has no molecular contrast, independent of its context.
    r=bank[ids[0]];a=tokenize(r['reference'][:51])[None,:].astype('int32');p=np.array([25],dtype='int32');h=np.array(encoder.encode(mx.array(a),mx.array(a),mx.array(p)))[:,128:];m=np.array(encoder.scan(mx.array(a),mx.array(a),mx.array(p)));assert np.max(abs(h))==0 and np.max(abs(m[:,155:310]))==0 and np.max(abs(patterns(a,a)[:,84:]))==0
    rng=np.random.default_rng(7701);c=rng.normal(0,.3,(14,560)).astype('float32');v=np.abs(rng.normal(0,.25,(14,623))).astype('float32');parameters=initialize(560,623,42);parameters['effect_w']=np.full((64,2),.1,dtype='float32');u=np.arange(14,dtype='int32')%2;target=np.linspace(.02,.8,14,dtype='float32');w=np.ones(14,dtype='float32');pm={k:mx.array(a) for k,a in parameters.items()};cpu=forward(parameters,c,v,backend='numpy');gpu=np.array(forward(pm,mx.array(c),mx.array(v)));error=float(np.max(abs(cpu-gpu)));assert error<2e-5
    assert np.max(abs(forward(parameters,c,np.zeros_like(v),backend='numpy')))==0
    masks=[mx.ones((14,128)),mx.ones((14,64)),mx.ones((14,128)),mx.ones((14,64))]
    def objective(pp):return loss(pp,mx.array(c),mx.array(v),mx.array(u),mx.array(target),mx.array(w),masks,6,2,2,.15)
    value,g=mx.value_and_grad(objective)(pm);direction=rng.normal(size=parameters['effect_w'].shape).astype('float32');direction/=np.linalg.norm(direction);eps=.001;plus=dict(pm);minus=dict(pm);plus['effect_w']=pm['effect_w']+eps*mx.array(direction);minus['effect_w']=pm['effect_w']-eps*mx.array(direction);numeric=(float(objective(plus))-float(objective(minus)))/(2*eps);analytic=float(np.sum(np.array(g['effect_w'])*direction));relative=abs(numeric-analytic)/max(abs(numeric),abs(analytic),1e-4);assert relative<.03,(numeric,analytic,relative)
    write_json(OUT/'feature_and_gradient_audit.json',{'passed':True,'independent_NumPy_feature_checks':checks,'forward_NumPy_MLX_max_error':error,'quantitative_loss_numeric_gradient':numeric,'quantitative_loss_backprop_gradient':analytic,'relative_gradient_error':relative,'null_contrast_features_and_magnitude_exactly_zero':True,'scope':'Source reconstruction and numerical implementation evidence, not biological predictive validity or clinical calibration.'});print('V7 INDEPENDENT FEATURES AND MAGNITUDE BACKPROP AUDIT PASSED',len(checks),relative,flush=True)
if __name__=='__main__':main()
