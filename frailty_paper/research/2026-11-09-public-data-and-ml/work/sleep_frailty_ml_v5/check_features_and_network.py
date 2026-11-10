"""Independent arithmetic, literal sequence-count and gradient checks."""
import sys,argparse,gzip,collections
from common import *
sys.path.insert(0,str(OLD));import model
from prepare_features import selected_motifs
import motif_model
from network import initialize,forward,branches,mx
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--partial-local-only',action='store_true');a=parser.parse_args();x=np.load(DATA/'fixed_feature_bank.npy',mmap_mode='r')
    with np.load(DATA/'local151_sequences.npz') as z:r=z['ref'][:4].astype('int32');s=z['alt'][:4].astype('int32');pos=z['position'][:4]
    groups,labels=selected_motifs();weights=model.read_weights();h=np.asarray(model.strand_features(r,s,pos,weights));motif=motif_model.scan_numpy(r,s,pos,groups)
    assert np.max(np.abs(h-x[:4,:512]))<5e-5;assert np.max(np.abs(motif-x[:4,1192:]))<5e-5
    def literal(a):
        text=np.frombuffer(b'ACGT',dtype='uint8')[a-7].tobytes().decode();out=[]
        for k in [1,2,3]:
            from itertools import product
            count=collections.Counter(text[i:i+k] for i in range(len(text)-k+1));out.extend(count[''.join(t)] for t in product('ACGT',repeat=k))
        return np.array(out,dtype='float64')
    den=np.concatenate([np.full(4**k,151-k+1) for k in [1,2,3]]);patterns=[]
    for aa,bb in zip(r,s):
        def one(a,b):
            p,q=literal(a),literal(b);return np.r_[.5*(p+q)/den,np.abs(p-q)]
        patterns.append(.5*(one(aa,bb)+one(17-aa[::-1],17-bb[::-1])))
    assert np.max(np.abs(np.stack(patterns)-x[:4,1024:1192]))<1e-6
    replay={'local_encoder_max_error':float(np.max(np.abs(h-x[:4,:512]))),'local_motif_max_error':float(np.max(np.abs(motif-x[:4,1192:]))),'literal_kmer_max_error':float(np.max(np.abs(np.stack(patterns)-x[:4,1024:1192])))}
    slices={'local_motif':[1192,x.shape[1]]};cols=np.arange(x.shape[1],dtype='int32');c,v=branches(cols,slices,3);rng=np.random.default_rng(9907);inp=rng.normal(0,.2,(8,x.shape[1]+3)).astype('float32');study=np.arange(8,dtype='int32')%3;y=np.array([0,1,0,1,0,0,1,0],dtype='float32');gradient=[]
    for architecture in ['mlp','contrast']:
        p=initialize(inp.shape[1],architecture,c,v,3,.1,942);pm={k:mx.array(a) for k,a in p.items()};xx=mx.array(inp);ss=mx.array(study);cc=mx.array(c);vv=mx.array(v)
        cpu=forward(p,inp,study,architecture,c,v,backend='numpy');gpu=forward(pm,xx,ss,architecture,cc,vv);error=max(float(np.max(np.abs(a-np.array(b)))) for a,b in zip(cpu,gpu));assert error<2e-5
        def objective(pp):
            z=forward(pp,xx,ss,architecture,cc,vv)[0];return mx.mean(mx.logaddexp(0,z)-mx.array(y)*z)+.2*mx.logaddexp(0,-(z[1]-z[0]))
        value,g=mx.value_and_grad(objective)(pm);key='class_w';direction=rng.normal(size=p[key].shape).astype('float32');direction/=np.linalg.norm(direction);eps=.002;plus=dict(pm);minus=dict(pm);plus[key]=pm[key]+eps*mx.array(direction);minus[key]=pm[key]-eps*mx.array(direction);numeric=(float(objective(plus))-float(objective(minus)))/(2*eps);analytic=float(np.sum(np.array(g[key])*direction));relative=abs(numeric-analytic)/max(abs(numeric),abs(analytic),1e-4);assert relative<.02,(architecture,numeric,analytic,relative)
        gradient.append({'architecture':architecture,'forward_max_error':error,'directional_numeric_gradient':numeric,'directional_backprop_gradient':analytic,'relative_difference':relative})
    if not a.partial_local_only:
        receipt=json.loads((OUT/'fixed_features_receipt.json').read_text());assert sha(receipt['feature_path'])==receipt['feature_sha256']
        with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
        sample=next(row for row in bank if row['sequence_source']=='kircher');i=int(sample['bank_id']);lut=np.zeros(256,dtype='uint8');lut[[ord(k) for k in 'ACGT']]=[7,8,9,10];rr=lut[np.frombuffer(sample['reference'].encode(),dtype='uint8')][None,:];aa=lut[np.frombuffer(sample['alternate'].encode(),dtype='uint8')][None,:];pp=np.array([int(sample['changed_position'])]);expected=np.asarray(model.strand_features(rr,aa,pp,weights));error=float(np.max(np.abs(expected-x[i:i+1,512:1024])));assert error<5e-5;replay['kircher_genomic_full_context_encoder_max_error']=error
        for start in range(0,len(x),2048):assert np.isfinite(x[start:start+2048]).all()
    report={'passed':True,'partial_local_only':a.partial_local_only,'fixed_feature_checks':replay,'independent_network_arithmetic_and_finite_difference_gradients':gradient,'scope':'Local NumPy versus MLX research-port checks and independent literal sequence counting; no claim of official external-runtime parity or biological accuracy.'}
    write_json(OUT/('partial_feature_gradient_check.json' if a.partial_local_only else 'full_feature_gradient_check.json'),report);print(json.dumps(report),flush=True)
if __name__=='__main__':main()
