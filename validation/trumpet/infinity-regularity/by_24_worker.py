"""Default Bowen--York free data at matched counts and identical H observers.

Use the already baseline-validated modal production backend to solve the
unchanged vanilla TwoPunctures equation. This is not a seed-geometry match:
BY P=m*Gamma*v and the Lorentz-transformed spatial spin match isolated
seed charges, not binary horizon charges.
"""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from bowen_york import Backend
from physical import constraints,norms

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--hispid-result',type=Path,required=True);p.add_argument('--observer',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execution',choices=('kokkos','reference'),default='kokkos');a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False)
r=json.loads(a.hispid_result.read_text());c=r['config'];b=Backend(str(a.library.resolve()));cfg=b.config();cfg['integer'].update(zip(('npoints_A','npoints_B','npoints_phi'),c['n']));cfg['real']['Newton_tol']=c['tolerance'];cfg['integer']['Newton_maxit']=24;cfg['integer']['verbose']=1
for side,h in zip(('plus','minus'),c['hole']):
 cfg['real']['par_m_'+side]=h['mass'];vel=np.array(h['velocity']);G=1/np.sqrt(1-vel@vel);momentum=h['mass']*G*vel
 rest_spin=np.array(h['spin']);spin=G*rest_spin-G*G/(G+1)*(vel@rest_spin)*vel
 for k in range(3):cfg['real']['par_P_'+side+str(k+1)]=float(momentum[k]);cfg['real']['par_S_'+side+str(k+1)]=float(spin[k])
assert np.array_equal(c['hole'][0]['center'],[3,0,0]) and np.array_equal(c['hole'][1]['center'],[-3,0,0])
assert hashlib.sha256(a.observer.read_bytes()).hexdigest()==r['observer_artifact_sha256']
with np.load(a.observer) as z:x=z['xyz'].copy();step=z['step'].copy()
result=dict(scope='same_BY_equation_as_vanilla_not_equal_HiSpID_free_metric',config=cfg,resolution=c['n'],
            hispid_result_sha256=hashlib.sha256(a.hispid_result.read_bytes()).hexdigest(),
            observer_sha256=r['observer_artifact_sha256'],library_sha256=b.library_sha256(),
            matching='bare masses, separation, Lorentz-transformed spin vectors, P=m*Gamma*v, collocation counts, fixed Cartesian H observer and FD steps; native coordinate maps and internal stopping weights differ')
(a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
start=time.monotonic()
with b.create(cfg) as s:
 result['diagnostics']=s.solve(linear_rtol=.001,krylov='gmres',preconditioner='modal',max_krylov=512,restart=64,execution=a.execution,memory_limit_mib=65536)
 result['solve_seconds']=time.monotonic()-start;result['work_statistics']=s.work_statistics()
 (a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 # Retain physical observers even when native roundoff prevents the declared
 # tolerance. Such a run remains explicitly unqualified.
 result['tolerance_gate_passed']=result['diagnostics']['converged']
 np.save(a.output/'unknowns.npy',s.unknowns())
 physical=constraints(s.sample,x,step)
 np.savez_compressed(a.output/'physical.npz',xyz=x,step=step,**physical)
 result['physical']={name:norms(physical,mask) for name,mask in [('near',np.arange(len(x))<48),('bulk',(np.arange(len(x))>=48)&(np.arange(len(x))<66))]}
result['physical_sha256']=hashlib.sha256((a.output/'physical.npz').read_bytes()).hexdigest();result['total_seconds']=time.monotonic()-start
(a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['physical']))
