"""Read-only current order-one boundary residuals on retained C1 states."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--executable',type=Path,required=True);p.add_argument('--runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False);rows=[]
for n in (160,192,224):
 binary=a.output/f'weights{n}.bin';q=subprocess.run([str(a.executable),str(n),str(binary)],check=True,capture_output=True,text=True)
 w=np.fromfile(binary).reshape(5,2,2*n*n).astype(np.longdouble);binary.unlink()
 root=a.runs/f'moderate{n}';r=json.loads((root/'result.json').read_text());assert hashlib.sha256((root/'solve.npz').read_bytes()).hexdigest()==r['solve_artifact_sha256']
 u=np.load(root/'solve.npz')['unknowns'].reshape(16,2*n*n,4).astype(np.longdouble);errors=[]
 for k in range(16):
  m=k if k<=8 else k-8
  if m<=4:errors.extend(np.einsum('sp,pv->sv',w[m],u[k]).ravel())
 error=float(max(abs(x) for x in errors));assert error<1e-12,error
 rows.append(dict(n=n,current_boundary_residual_linf=error,weight_comparison=q.stderr.strip(),solve_sha256=r['solve_artifact_sha256']))
result=dict(passed=True,scope=__doc__,new_elliptic_solve=False,physical_acceptance_inherited=False,records=rows)
(a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
