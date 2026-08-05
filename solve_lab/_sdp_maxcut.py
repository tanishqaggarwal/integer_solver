"""Order-1 (Shor/MAX-CUT) SDP relaxation + Goemans-Williamson rounding on the reduced
boolean problem. Ground truth: exact fail-count f(b) over all 2^10 combos of the 10
residual-feeding selectors (scan10.json). Fit quadratic model, relax, GW-round many seeds."""
import json, numpy as np, itertools, cvxpy as cp
SEL10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
raw=json.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/scan10.json'))
# raw: list of (nfail, combo(list10), failidxs)
data={tuple(c):n for n,c,f in raw}
n=10
# ground truth optimum
gt=min(data.values()); argmins=[c for c,v in data.items() if v==gt]
print(f'ground-truth min fails over 2^{n}={len(data)} combos: {gt}  at {len(argmins)} pts e.g. {argmins[0]}')
# Fit quadratic model f(b) ~ c0 + sum g_i b_i + sum_{i<j} Q_ij b_i b_j  (b in {0,1})
rows=[]; y=[]
idx=[(i,j) for i in range(n) for j in range(i+1,n)]
for c,v in data.items():
    r=[1.0]+[c[i] for i in range(n)]+[c[i]*c[j] for (i,j) in idx]
    rows.append(r); y.append(v)
A=np.array(rows); y=np.array(y,float)
coef,res,rk,sv=np.linalg.lstsq(A,y,rcond=None)
pred=A@coef
print(f'quadratic-model fit: R^2={1-np.sum((pred-y)**2)/np.sum((y-y.mean())**2):.4f}, maxabs resid={np.max(np.abs(pred-y)):.2f}')
c0=coef[0]; g=coef[1:1+n]; Qd={}
for k,(i,j) in enumerate(idx): Qd[(i,j)]=coef[1+n+k]
# Build symmetric Q for x in {-1,1} via b=(1+x)/2 -> map quadratic in b to quadratic in x.
# f(b)=c0+ g.b + sum Qij bi bj ;  bi=(1+xi)/2
# We minimize; encode as MAX-CUT-style: minimize x^T W x + h^T x over x in {-1,1}
# Expand:
W=np.zeros((n,n)); h=np.zeros(n); const=c0
for i in range(n):
    const+=g[i]*0.5; h[i]+=g[i]*0.5
for (i,j),q in Qd.items():
    # bi bj = (1+xi)(1+xj)/4 = 1/4(1 + xi + xj + xi xj)
    const+=q*0.25; h[i]+=q*0.25; h[j]+=q*0.25; W[i,j]+=q*0.25; W[j,i]+=q*0.25
# homogenize with an extra spin x0 (=+1): minimize [x;x0]^T M [x;x0]
M=np.zeros((n+1,n+1)); M[:n,:n]=W
for i in range(n): M[i,n]=h[i]/2; M[n,i]=h[i]/2
# SDP relaxation: X psd, diag(X)=1, minimize <M,X>
X=cp.Variable((n+1,n+1),symmetric=True)
cons=[X>>0]+[X[i,i]==1 for i in range(n+1)]
prob=cp.Problem(cp.Minimize(cp.trace(M@X)+const),cons)
prob.solve(solver=cp.SCS,verbose=False)
print(f'SDP relaxation lower bound (model units): {prob.value:.3f}  (status {prob.status})')
Xv=X.value
# GW rounding many hyperplanes
np.random.seed(0)
w,V=np.linalg.eigh((Xv+Xv.T)/2)
w=np.clip(w,0,None)
L=(V*np.sqrt(w))  # rows are vectors; L@r gives embedding
best=(1e9,None)
for t in range(20000):
    r=np.random.randn(n+1); s=np.sign(L@r)
    s=s*s[n]  # fix anchor spin to +1
    b=((1+s[:n])/2).astype(int)
    fval=data.get(tuple(b.tolist()),None)
    if fval is not None and fval<best[0]:
        best=(fval,tuple(b.tolist()))
print(f'GW rounding (20000 hyperplanes): best TRUE fails found = {best[0]}  at {best[1]}')
on=[SEL10[i] for i in range(n) if best[1][i]]
print(f'  -> selectors ON: {on}')
print(f'CONCLUSION: SDP+GW best={best[0]} vs ground-truth={gt}; improvement over 11? {best[0]<11}')
