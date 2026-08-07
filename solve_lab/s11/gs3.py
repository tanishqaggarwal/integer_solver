import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2, gs
P=L.P
MATCH=[(688,30213),(1618,22162),(26719,14515),(26721,19750),(26733,5096),(28438,21589)]
WATCH=[688,1618,26719,26721,26723,26733,28438,32342,36185,40608]
theta={c:0 for _,c in MATCH}
def resid(th, atoms):
    v=gs.state(th); return v,[fw.evalpoly(L.polys[a],v)%P for a in atoms]
t0=time.time()
for it in range(14):
    v,r=resid(theta, WATCH)
    nz=[WATCH[i] for i,x in enumerate(r) if x!=0]
    print(f"it{it}: nonzero mod p: {nz}  ({time.time()-t0:.0f}s)")
    if not nz: break
    for a,c in MATCH:
        v,r1=resid(theta,[a])
        if r1[0]==0: continue
        th=dict(theta); th[c]=theta[c]+1
        _,r2=resid(th,[a])
        s=(r2[0]-r1[0])%P
        if s==0: continue
        theta[c]=(theta[c] + (-r1[0])*pow(s,-1,P)) % P
v,r=resid(theta,WATCH)
print("final mod-p residuals zero:", all(x==0 for x in r))
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(k):val for k,val in theta.items()}, open('theta_modp.json','w'))
