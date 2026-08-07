import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine
from zsolve import solve_int
CTRL=json.load(open('controls.json'))
def resid(theta, BAD):
    v=engine.apply_theta(theta)
    return [fw.evalpoly(L.polys[a],v) for a in BAD], v
BAD=[26719,26721,26723,26733,28438,32342,36185]
theta={c:0 for c in CTRL}
t0=time.time()
r0,v0=resid(theta,BAD)
print("base residual magnitudes:", [len(str(abs(x))) for x in r0])
J=[[0]*len(CTRL) for _ in BAD]
lin=True
for j,c in enumerate(CTRL):
    th=dict(theta); th[c]=1
    r1,_=resid(th,BAD)
    th[c]=2
    r2,_=resid(th,BAD)
    for i in range(len(BAD)):
        d1=r1[i]-r0[i]; d2=r2[i]-r0[i]
        J[i][j]=d1
        if d2!=2*d1: lin=False; print(f"  NONLINEAR: check a{BAD[i]} wrt x{c}")
print(f"jacobian built ({time.time()-t0:.0f}s) linear={lin}")
for i,a in enumerate(BAD):
    print(f"  a{a}: nz cols {[(CTRL[j], len(str(abs(J[i][j])))) for j in range(len(CTRL)) if J[i][j]]}")
x=solve_int(J,[-t for t in r0])
print("integer solution:", "FOUND" if x else "NONE")
if x:
    th={CTRL[j]: x[j] for j in range(len(CTRL))}
    print({k:str(val)[:24] for k,val in th.items() if val})
    json.dump({str(k):val for k,val in th.items()}, open('theta.json','w'))
    r,v=resid(th,BAD)
    print("residuals after:", [str(t)[:12] for t in r])
    b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
    print("bad:", b)
    json.dump({str(i):v[i] for i in range(L.NVARS)}, open('newton2.json','w'))
