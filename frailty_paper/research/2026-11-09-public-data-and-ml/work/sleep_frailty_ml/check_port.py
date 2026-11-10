from pathlib import Path
import sys,json,time
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np
import mlx.core as mx
import model as reference
import model_mlx as model
from autograd import grad

OUT=Path(__file__).resolve().parents[2]/'outputs/sleep_frailty_ai_v1'
OUT.mkdir(exist_ok=True)
npw=reference.read_weights();w={k:mx.array(v) for k,v in npw.items()}
r=np.random.default_rng(25)
tokens=r.integers(7,11,(2,32));expected=reference.encode(tokens,npw);actual=np.array(model.encode(mx.array(tokens),w))
error=float(np.max(np.abs(expected-actual)));assert error<2e-4,error
u=r.normal(size=(2,3,32)).astype('float32');k=r.normal(size=(3,32)).astype('float32');d=r.normal(size=3).astype('float32')
convolution_error=float(np.max(np.abs(np.array(model.fftconv(mx.array(u),mx.array(k),mx.array(d)))-reference.direct_conv(u,k,d))))
assert convolution_error<1e-4
conv_loss=lambda kk:reference.np.sum(reference.fftconv(u,kk,d)**2)
expected_k_grad=grad(conv_loss)(k.astype(float))
actual_k_grad=np.array(mx.grad(lambda kk:mx.sum(model.fftconv(mx.array(u),kk,mx.array(d))**2))(mx.array(k)))
conv_gradient_error=float(np.max(np.abs(expected_k_grad-actual_k_grad)))
assert np.allclose(actual_k_grad,expected_k_grad,atol=1e-3,rtol=1e-4),conv_gradient_error
# Check gradients spanning FFT/filter, linear, normalization and feed-forward paths.
checks=[]
for key,index in [('layers.1.mixer.in_proj.weight',(0,0)),('layers.1.mixer.filter_fn.implicit_filter.0.weight',(0,0)),('layers.1.mlp.fc2.bias',(0,)),('ln_f.weight',(0,))]:
    v=npw[key].astype(float)
    fn=lambda a:reference.encode(tokens,{**npw,key:a})[0,0,0]
    ag=float(grad(fn)(v)[index]);eps=1e-4
    plus=v.copy();minus=v.copy();plus[index]+=eps;minus[index]-=eps
    fd=float((fn(plus)-fn(minus))/(2*eps))
    mf=lambda a:model.encode(mx.array(tokens),{**w,key:a})[0,0,0]
    mg=float(mx.grad(mf)(w[key])[index])
    assert abs(ag-fd)<2e-4+abs(fd)*.005,(key,ag,fd)
    assert abs(mg-ag)<2e-4+abs(ag)*.005,(key,mg,ag)
    checks.append({'parameter':key,'autograd':ag,'finite_difference':fd,'mlx':mg})
L=mx.arange(8,dtype=mx.float32)
real_bad=np.array(mx.grad(lambda a:mx.sum(mx.fft.irfft(mx.fft.rfft(a))**2))(L))
complex_good=np.array(mx.grad(lambda a:mx.sum(mx.real(mx.fft.ifft(mx.fft.fft(a)))**2))(L))
assert np.max(np.abs(complex_good-2*np.arange(8)))<1e-4
test=r.integers(7,11,(16,270));alternate=test.copy();alternate[:,134]=7+(alternate[:,134]-6)%4
t=time.time();features=model.strand_features(mx.array(test),mx.array(alternate),mx.full((16,),134),w);mx.eval(features);seconds=time.time()-t
# Swapping alleles and reversing-complementing the pair must preserve detection features.
swap=model.strand_features(mx.array(alternate),mx.array(test),mx.full((16,),134),w)
reverse=model.strand_features(mx.array(17-test[:,::-1]),mx.array(17-alternate[:,::-1]),mx.full((16,),135),w)
swap_error=float(mx.max(mx.abs(features-swap)));rc_error=float(mx.max(mx.abs(features-reverse)))
assert swap_error<1e-6 and rc_error<1e-6
receipt={'encoder_forward_max_error':error,'convolution_direct_max_error':convolution_error,'convolution_parameter_gradient_max_error':conv_gradient_error,'parameter_gradient_checks':checks,'complex_fft_identity_gradient_max_error':float(np.max(np.abs(complex_good-2*np.arange(8)))),'unused_real_fft_identity_gradient_max_error':float(np.max(np.abs(real_bad-2*np.arange(8)))),'allele_swap_max_error':swap_error,'reverse_complement_max_error':rc_error,'batch16_pair_feature_seconds':seconds,'feature_dim':features.shape[1],'device':str(mx.default_device()),'official_pytorch_executable_parity':'Not tested; source-formula NumPy/Autograd port and independent convolution/gradient checks used instead'}
(OUT/'port_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
