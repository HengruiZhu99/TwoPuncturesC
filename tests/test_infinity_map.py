"""Public sampler independently checked against a manufactured mapped field."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'python'))
from hispid import Backend
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();b=Backend(str(a.library.resolve()));assert any(k in b.parameterization() for k in ('inf2_map_v6','poly2_map_v7','beta14_map_v8','focus_map_v9','focus2_map_v10'))
c=b.config();c.far_radius=0;c.n[:]=[12,24,16];c.hole[0].center[:]=[3,0,0];c.hole[1].center[:]=[-3,0,0]
na,nb,np_=c.n;maps=b.parameterization_maps();l=maps['radial_stretch'];kappa=maps['angular_stretch']
s=(1-np.cos(np.pi*(np.arange(na)+.5)/na))/2;z=-np.cos(np.pi*(np.arange(nb)+.5)/nb)
u=np.zeros((np_,nb,na,4));mode=3
u[mode,:,:,0]=1e-4*(1+.2*s[None,:]+.07*z[:,None]+.03*s[None,:]*z[:,None]+.02*s[None,:]**2)
# Analytic prolate coordinates, independent of native inverse-map code.
t=np.array([.002,.07,.3,.8,.98]);eta=np.array([-.91,.67,.3,-.4,.1]);phi=np.array([.1,.5,1.1,2.,3.])
x=np.stack([3*(1+t)/(1-t)*eta,6*np.sqrt(t)/(1-t)*np.sqrt(1-eta*eta)*np.cos(phi),6*np.sqrt(t)/(1-t)*np.sqrt(1-eta*eta)*np.sin(phi)],axis=1)
def beta_inverse(tt):
 lo=np.zeros_like(tt.real);hi=np.ones_like(tt.real);bb=(2-l)/5;cc=1-l-bb
 def forward(ss):
  qq=1-ss;return l*ss+bb*ss**5+cc*ss**15*(136*qq*qq+17*ss*qq+ss*ss)
 for _ in range(60):
  mid=(lo+hi)/2;lower=forward(mid)<tt.real;lo=np.where(lower,mid,lo);hi=np.where(lower,hi,mid)
 ss=((lo+hi)/2).astype(tt.dtype)
 for _ in range(8):
  qq=1-ss;derivative=l+5*bb*ss**4+2040*cc*ss**14*qq*qq
  ss=ss-(forward(ss)-tt)/derivative
 return ss
q=l*(1-t)
if maps.get('radial_kind')=='endpoint_beta14':s_eval=beta_inverse(t)
elif maps.get('radial_kind')=='quadratic':s_eval=2*t/(l+np.sqrt(l*l+4*(1-l)*t))
elif 'inf2_map_v6' in b.parameterization():s_eval=2*t/(2*t+q+np.sqrt(q*(4*t+q)))
else:s_eval=t/(l+(1-l)*t)
z_eval=np.arctanh(eta*np.tanh(kappa))/kappa
P=1e-4*(1+.2*s_eval+.07*z_eval+.03*s_eval*z_eval+.02*s_eval*s_eval)
expected=-2*(1-t)*(t*(1-eta*eta))**(mode/2)*P*np.sqrt(2/np_)*np.cos(mode*phi)
with b.create_sampler(c) as sol:
 sol.set_unknowns(u.ravel());sampled=sol.sample_with_derivatives(x);got=sampled['correction'][:,0]
err=float(np.max(abs(got-expected)));assert err<2e-16,(got,expected)
# Independent Cartesian complex-step oracle. The rest seeds are two
# Schwarzschild holes, so gamma=(1+.25/rp+.25/rm+u)^4 delta exactly.
def analytic_gamma(points):
 X,Y,Z=points.T;rp=np.sqrt((X-3)**2+Y*Y+Z*Z);rm=np.sqrt((X+3)**2+Y*Y+Z*Z)
 D=(rp+rm+6)/2;tt=1-6/D;ee=(rm-rp)/6
 if maps.get('radial_kind')=='endpoint_beta14':ss=beta_inverse(tt)
 elif maps.get('radial_kind')=='quadratic':ss=2*tt/(l+np.sqrt(l*l+4*(1-l)*tt))
 elif 'inf2_map_v6' in b.parameterization():
  qq=l*(1-tt);ss=2*tt/(2*tt+qq+np.sqrt(qq*(4*tt+qq)))
 else:ss=tt/(l+(1-l)*tt)
 zz=np.arctanh(ee*np.tanh(kappa))/kappa
 PP=1e-4*(1+.2*ss+.07*zz+.03*ss*zz+.02*ss*ss)
 correction=-12*np.sqrt(2/np_)*(Y**3-3*Y*Z*Z)*PP/D**4
 psi=1+.25/rp+.25/rm+correction
 return psi[:,None,None]**4*np.eye(3)
expected_gradient=np.zeros((len(x),3,3,3));delta=1e-20
for d in range(3):
 xx=x.astype(complex);xx[:,d]+=1j*delta
 expected_gradient[:,d]=analytic_gamma(xx).imag/delta
gradient_error=float(np.max(abs(sampled['dgamma']-expected_gradient)/(1+abs(expected_gradient))))
assert gradient_error<2e-11,gradient_error
a.output.write_text(json.dumps(dict(passed=True,absolute_error=err,gradient_scaled_linf=gradient_error,points=x.tolist(),library_sha256=b.library_sha256(),parameterization=b.parameterization()),indent=2)+'\n');print(err)
