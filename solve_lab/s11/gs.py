"""Mod-p Gauss-Seidel: match each bad check to a control it is LINEAR in, solve, iterate."""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2
P=L.P
NAT={u:len(L.var_atoms[u]) for u in range(L.NVARS)}
B4=(542,47,438,91)

def state(theta, derive=False):
    return engine2.close(B4, theta, derive=derive)

def probe(theta, C, BAD, derive=False):
    """return base residuals and slope[b][c] mod p (None if nonlinear)"""
    v=state(theta, derive)
    base=[fw.evalpoly(L.polys[a],v)%P for a in BAD]
    slope={}
    for c in C:
        th=dict(theta)
        th[c]=theta.get(c,0)+1; v1=state(th,derive); r1=[fw.evalpoly(L.polys[a],v1)%P for a in BAD]
        th[c]=theta.get(c,0)+2; v2=state(th,derive); r2=[fw.evalpoly(L.polys[a],v2)%P for a in BAD]
        for i,a in enumerate(BAD):
            d1=(r1[i]-base[i])%P; d2=(r2[i]-r1[i])%P
            if d1==0 and d2==0: continue
            slope[(a,c)]= d1 if d1==d2 else None
    return v, base, slope

if __name__=='__main__':
    t0=time.time()
    v=state({})
    b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"no-derive 4bit: bad={len(b)} failing={len(f)} score={L.NEQ-len(f)} ({time.time()-t0:.0f}s)")
    print("bad:", b)
    modp=[a for a in b if fw.evalpoly(L.polys[a],v)%P!=0]
    print("bad mod p:", modp)
