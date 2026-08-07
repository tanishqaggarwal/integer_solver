"""Cascade closer: zero every nonzero check atom by assigning a DIRECT free input of it."""
import ev, fast, json, os, sys, time, math, pickle, random
from fast import St, csup, inv, chk
HERE=os.path.dirname(os.path.abspath(__file__))
FREE=set(ev.F['free0'])
AVARS=ev.atom_vars

def coef_and_rest(st,a,X):
    """value of atom a as a function of X: probe at X, X+1, X+2 (deg<=2 in X)."""
    base=st.fv.get(X,0)
    t=st.clone()
    y0=st.av[a]
    t.set_free({X:base+1}); y1=t.av[a]
    t.set_free({X:base+2}); y2=t.av[a]
    d1=y1-y0; d2=y2-2*y1+y0
    return y0,d1,d2,base

def roots(st,a,X):
    y0,d1,d2,base=coef_and_rest(st,a,X)
    if d2==0:
        if d1==0: return []
        if y0%d1: return []
        return [base-y0//d1]
    # y(t)=y0 + d1*t + d2*t(t-1)/2
    aa=d2; bb=2*d1-d2; cc=2*y0
    disc=bb*bb-4*aa*cc
    if disc<0: return []
    r=math.isqrt(disc)
    if r*r!=disc: return []
    out=[]
    for s in (-bb+r,-bb-r):
        if s%(2*aa)==0: out.append(base+s//(2*aa))
    return out

def close(st, maxsteps=400, verbose=False, forbid=frozenset(), rng=None):
    st=st.clone()
    for step in range(maxsteps):
        nz=st.nz()
        if not nz: return st,True
        best=None
        for a in nz:
            for X in AVARS[a]:
                if X not in FREE or X in forbid: continue
                for rt in roots(st,a,X):
                    if st.fv.get(X,0)==rt: continue
                    g=st.clone().set_free({X:rt})
                    if g.av[a]!=0: continue
                    k=(len(g.nz()),len(g.fails))
                    if best is None or k<best[0]: best=(k,g,a,X,rt)
        if best is None: return st,False
        k,g,a,X,rt=best
        if verbose: print('  a%d <- x_%d  nz=%d fails=%d'%(a,X,k[0],k[1]),flush=True)
        st=g
    return st,False

if __name__=='__main__':
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    ub=int(sys.argv[1]); wb=int(sys.argv[2])
    fv={5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,ub:1,wb:1}
    st=St(fv)
    print('start score',st.score(),'nz',sorted(st.nz()))
    t0=time.time()
    out,ok=close(st,verbose=True)
    print('closed=%s score=%d nz=%s  %.1fs'%(ok,out.score(),sorted(out.nz()),time.time()-t0))
    if out.score()>=39005:
        p=os.path.join(HERE,'C_%d_%d_%d.json'%(out.score(),ub,wb))
        json.dump({('x_%d'%i):out.v[i] for i in range(38748) if out.v[i]!=0},open(p,'w'))
        json.dump({str(a):b for a,b in out.fv.items()},open(p.replace('.json','_fv.json'),'w'))
        print('wrote',p)
