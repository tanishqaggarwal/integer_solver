"""Parametrised structural-closure engine."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P=L.P
NAT={u:len(L.var_atoms[u]) for u in range(L.NVARS)}
C0=L.polys[688][()]; MM=8863713; G0=(-C0*pow(MM,-1,P))%P; C0B=L.polys[1618][()]
DERIVED={14853,16742,14515,21589,30213,22820,7497,22162,14393,11436}
def free_cands(a, locked):
    out=[]
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked: continue
        if any(mm.count(u)>1 for mm in L.polys[a]): continue
        out.append(u)
    out.sort(key=lambda u:(NAT[u],u)); return out
def close(BITS, theta, pin_rounds=25, derive=True):
    v=[0]*L.NVARS
    for b in BITS: v[b]=1
    for k,x in theta.items(): v[k]=x
    fw.forward(v)
    if derive:
        v[14853]=v[12186]; v[16742]=v[24908]
        fw.forward(v)
        v[14515]+= v[14853]-v[1308]
        v[21589]+= v[16742]-v[19083]
        fw.forward(v)
        v[30213]=G0; v[22820]=0; v[7497]=(C0+MM*G0)//P
        v[22162]=-C0B; v[14393]=0; v[11436]=0
        fw.forward(v)
    locked=set(BITS)|(DERIVED if derive else set())|set(theta)
    for it in range(pin_rounds):
        b=fw.bad_checks(v)
        if not b: break
        prog=False
        for a in b:
            if fw.evalpoly(L.polys[a],v)==0: continue
            for t in free_cands(a, locked):
                x=fw.solve_lin(a,t,v)
                if x is not None and x!=v[t]:
                    old=v[t]; v[t]=x; fw.forward(v)
                    if fw.evalpoly(L.polys[a],v)==0: prog=True; break
                    v[t]=old; fw.forward(v)
        if not prog: break
    return v
