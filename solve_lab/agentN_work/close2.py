"""Bottom-up cascade closer: smallest-support atom first, freeze assigned variables."""
import ev, fast, json, os, sys, time, math
from fast import St, csup, inv
from close import roots
HERE=os.path.dirname(os.path.abspath(__file__))
FREE=set(ev.F['free0']); AVARS=ev.atom_vars
def popc(x): return bin(x).count('1')
SUPN={a:popc(s) for a,s in csup.items()}

def close(st, maxsteps=600, verbose=False, frozen=None):
    st=st.clone()
    frozen=set(frozen or ())
    for step in range(maxsteps):
        nz=st.nz()
        if not nz: return st,True,frozen
        nz.sort(key=lambda a: SUPN[a])
        moved=False
        for a in nz:
            cands=[X for X in AVARS[a] if X in FREE and X not in frozen]
            cands.sort(key=lambda X: popc(csup.get(a,0)))
            best=None
            for X in cands:
                for rt in roots(st,a,X):
                    if st.fv.get(X,0)==rt: continue
                    g=st.clone().set_free({X:rt})
                    if g.av[a]!=0: continue
                    k=(len(g.nz()),len(g.fails))
                    if best is None or k<best[0]: best=(k,g,X,rt)
            if best is None: continue
            k,g,X,rt=best
            st=g; frozen.add(X); moved=True
            if verbose: print('  step%d a%d(sup%d) <- x_%d : nz=%d fails=%d'%(step,a,SUPN[a],X,k[0],k[1]),flush=True)
            break
        if not moved: return st,False,frozen
    return st,False,frozen

if __name__=='__main__':
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    ub=int(sys.argv[1]); wb=int(sys.argv[2])
    fv={5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,ub:1,wb:1}
    st=St(fv)
    print('start score',st.score(),'nz',sorted(st.nz()))
    t0=time.time()
    # freeze the deliberately-set values
    out,ok,fr=close(st,verbose=True,frozen=set(fv))
    print('closed=%s score=%d nz=%s  %.1fs'%(ok,out.score(),sorted(out.nz()),time.time()-t0))
    p=os.path.join(HERE,'D_%d_%d_%d.json'%(out.score(),ub,wb))
    json.dump({('x_%d'%i):out.v[i] for i in range(38748) if out.v[i]!=0},open(p,'w'))
    json.dump({str(a):b for a,b in out.fv.items()},open(p.replace('.json','_fv.json'),'w'))
    print('wrote',p)
