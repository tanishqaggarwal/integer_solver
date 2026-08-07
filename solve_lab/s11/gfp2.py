import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, gs, gfp
P=L.P
theta={c:0 for c in gfp.CTRL}
v,r=gfp.resid(theta)
J=[[0]*len(gfp.CTRL) for _ in gfp.TARGETS]
for j,c in enumerate(gfp.CTRL):
    th=dict(theta); th[c]=1
    _,r1=gfp.resid(th)
    for i in range(len(gfp.TARGETS)): J[i][j]=(r1[i]-r[i])%P
print("controls:", gfp.CTRL)
for i,(nm,_) in enumerate(gfp.TARGETS):
    print(f"  {nm:8s}: affected by {[gfp.CTRL[j] for j in range(len(gfp.CTRL)) if J[i][j]]}")
