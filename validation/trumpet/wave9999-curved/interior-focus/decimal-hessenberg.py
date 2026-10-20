"""High-precision replay of only the recorded small Hessenberg least squares."""
from decimal import Decimal,localcontext
from pathlib import Path
import argparse,gzip,json,time
parser=argparse.ArgumentParser();parser.add_argument('input',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
with (gzip.open(args.input,'rt') if args.input.suffix=='.gz' else args.input.open()) as stream:r=json.loads(next(stream))
assert r['method']==2 and r['cycle_begin']==0 and 1<=r['columns']<=256
n=r['columns'];raw=r['hessenberg']
results=[]
for digits in (35,65):
 start=time.monotonic()
 with localcontext() as c:
  c.prec=digits
  h=[[Decimal(x) for x in row] for row in raw];g=[Decimal(0)]*(n+1);g[0]=Decimal(1);cs=[];sn=[]
  for k in range(n):
   for j in range(k):
    a=cs[j]*h[j][k]+sn[j]*h[j+1][k]
    h[j+1][k]=-sn[j]*h[j][k]+cs[j]*h[j+1][k];h[j][k]=a
   dd=(h[k][k]**2+h[k+1][k]**2).sqrt();cs.append(h[k][k]/dd);sn.append(h[k+1][k]/dd);h[k][k]=dd;h[k+1][k]=Decimal(0)
   g[k+1]=-sn[k]*g[k];g[k]=g[k]*cs[k]
  y=[Decimal(0)]*n
  for j in range(n-1,-1,-1):y[j]=(g[j]-sum((h[j][k]*y[k] for k in range(j+1,n)),Decimal(0)))/h[j][j]
  res=[(Decimal(1) if i==0 else Decimal(0))-sum((Decimal(raw[i][j])*y[j] for j in range(n)),Decimal(0)) for i in range(n+1)]
  exact=sum((v*v for v in res),Decimal(0)).sqrt()
  result=dict(decimal_digits=digits,recurrence_relative=str(abs(g[-1])),projected_actual_relative=str(exact),coefficient_max=str(max(map(abs,y))),seconds=time.monotonic()-start)
  results.append(result);print(result,flush=True)
args.output.write_text(json.dumps(dict(scope='Exact rounded recorded H only; no full operator/solution acceptance',results=results),indent=2)+'\n')
