import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2
from gfp import gauss_solve
P=L.P
BITS=(542,47,438,91)
TARGETS=[('x3719',lambda v:v[3719]),('x25118',lambda v:v[25118]),
         ('x25614',lambda v:v[25614]),('x34220',lambda v:v[34220]),
         ('n-gap',lambda v:v[12186]-v[1308]),('m-gap',lambda v:v[24908]-v[19083])]
CTRL=[14515,19750,5096,21589,33708,31339,28486]
def resid(th):
    v=engine2.close(BITS, th, derive=False)
    return v,[f(v)%P for _,f in TARGETS]
theta={c:0 for c in CTRL}
t0=time.time()
for it in range(15):
    v,r=resid(theta)
    nz=[TARGETS[i][0] for i,x in enumerate(r) if x]
    print(f"it{it}: nonzero={nz} ({time.time()-t0:.0f}s)")
    if not nz: break
    J=[[0]*len(CTRL) for _ in TARGETS]
    for j,c in enumerate(CTRL):
        th=dict(theta); th[c]=theta[c]+1
        _,r1=resid(th)
        for i in range(len(TARGETS)): J[i][j]=(r1[i]-r[i])%P
    d=gauss_solve(J,[(-x)%P for x in r],P)
    if d is None:
        print("  inconsistent linearisation"); break
    for j,c in enumerate(CTRL): theta[c]=(theta[c]+d[j])%P
v,r=resid(theta)
print("ALL TARGETS ZERO:", not any(r))
if not any(r):
    json.dump({str(k):val for k,val in theta.items()}, open('theta_ok.json','w'))
    b=fw.bad_checks(v); print("bad checks:",len(b))
    print("mod-p nonzero:", [a for a in b if fw.evalpoly(L.polys[a],v)%P!=0])
