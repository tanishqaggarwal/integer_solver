import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2, gs
P=L.P
# triangular targets: (target VAR to zero mod p, control free input)
STEPS=[(3719,14515),(25118,19750),(25614,5096),(34220,21589),(None,None)]
theta={14515:0,19750:0,5096:0,21589:0,22162:0,30213:0}
def st(th): return gs.state(th)
def val(v,var): return v[var]%P
t0=time.time()
for rnd in range(6):
    for tgt,c in STEPS:
        if tgt is None: continue
        v=st(theta); r=val(v,tgt)
        if r==0: continue
        th=dict(theta); th[c]=theta[c]+1
        s=(val(st(th),tgt)-r)%P
        if s==0: print(f"  target x{tgt}: control x{c} has ZERO slope"); continue
        theta[c]=(theta[c]+(-r)*pow(s,-1,P))%P
    v=st(theta)
    print(f"rnd{rnd}: x3719={val(v,3719)!=0} x25118={val(v,25118)!=0} x25614={val(v,25614)!=0} x34220={val(v,34220)!=0} ({time.time()-t0:.0f}s)")
    if all(val(v,t)==0 for t,_ in STEPS if t): break
# now a688 / a1618
for a,c in [(688,30213),(1618,22162)]:
    v=st(theta); r=fw.evalpoly(L.polys[a],v)%P
    if r:
        th=dict(theta); th[c]=theta[c]+1
        s=(fw.evalpoly(L.polys[a],st(th))%P - r)%P
        if s: theta[c]=(theta[c]+(-r)*pow(s,-1,P))%P
v=st(theta)
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
print("mod-p nonzero:", [a for a in b if fw.evalpoly(L.polys[a],v)%P!=0])
json.dump({str(k):val_ for k,val_ in theta.items()}, open('theta_tri.json','w'))
