"""Independent scalar root/derivative and allele-symmetry checks before training."""
from measurement_common import *
from folded_effect_network import *
from scipy.optimize import brentq
from scipy.special import ndtr

def main():
    assert_prior_versions();delta=np.array([0.,1e-4,.02,.2,1.,-1.,2.,-3.,10.]);sigma=np.array([1.,.03,.1,.2,1.,1.,.05,.3,.01]);expected=np.array([brentq(lambda x:ndtr((x-d)/s)+ndtr((x+d)/s)-1.5,0,abs(d)+10*s,xtol=1e-14) for d,s in zip(delta,sigma)])
    d=torch.tensor(delta,dtype=torch.float64,requires_grad=True);s=torch.tensor(sigma,dtype=torch.float64,requires_grad=True);got=observed_magnitude_median(d,s);error=float(abs(got.detach().numpy()-expected).max());assert error<1e-8;got.sum().backward();checks=[]
    for i,(dd,ss) in enumerate(zip(delta,sigma)):
        epsilon=min(1e-5,ss*1e-4)
        def median(a,b):return brentq(lambda x:ndtr((x-a)/b)+ndtr((x+a)/b)-1.5,0,abs(a)+10*b,xtol=1e-14)
        gd=(median(dd+epsilon,ss)-median(dd-epsilon,ss))/(2*epsilon);gs=(median(dd,ss+epsilon)-median(dd,ss-epsilon))/(2*epsilon);difference=max(abs(gd-float(d.grad[i])),abs(gs-float(s.grad[i])));assert difference<1e-6,(i,difference);checks.append(difference)
    torch.manual_seed(42);model=FoldedEffectHead([.1,.2],[.1,.2],[0,1]).double().eval()
    with torch.no_grad():
        model.delta_output.weight.fill_(.03);model.mean_output.weight.fill_(.02);model.noise_output.weight.fill_(.01)
    generator=torch.Generator().manual_seed(971);mean=torch.randn(7,96,generator=generator,dtype=torch.float64);diff=torch.randn(7,96,generator=generator,dtype=torch.float64);physical=torch.randn(7,2,generator=generator,dtype=torch.float64);rr=torch.randn(7,3,generator=generator,dtype=torch.float64);aa=torch.randn(7,3,generator=generator,dtype=torch.float64);cell=torch.tensor([0,1,0,1,0,1,0]);flags=F.one_hot(cell,2).double();first=model(mean,diff,physical,rr,aa,cell,flags,flags);swapped=model(mean,-diff,physical,aa,rr,cell,flags,flags);zero=model(mean,torch.zeros_like(diff),physical,rr,rr,cell,flags,flags)
    assert torch.max(abs(first['signed_effect']+swapped['signed_effect']))<1e-10 and torch.max(abs(first['magnitude']-swapped['magnitude']))<1e-10;assert torch.max(abs(zero['signed_effect']))<1e-10
    write_json(DEST/'folded_median_pretraining_audit.json',{'passed':True,'independent_SciPy_median_max_error':error,'independent_implicit_derivative_max_error':max(checks),'allele_swap_signed_antisymmetry':True,'allele_swap_observed_magnitude_symmetry':True,'identical_DNA_zero_signed_effect':True,
              'interpretation':'Numerical and structural correctness only. Gaussian noise and supplied SE are working assumptions, not established biological truth or transferable confidence.',
              'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('folded_effect_network.py')]}});assert_prior_versions();print('FOLDED OBSERVATION MEDIAN AND ALLELE INVARIANTS PASSED',flush=True)

if __name__=='__main__':main()
