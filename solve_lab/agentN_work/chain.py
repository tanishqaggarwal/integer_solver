"""Enumerate the cascade's pin atoms: the ordered chain close2 uses to close a one-selector state."""
import ev, fast, json
from fast import St, csup
from close import roots
FREE=set(ev.F['free0']); AVARS=ev.atom_vars
def popc(x): return bin(x).count('1')
SUPN={a:popc(s) for a,s in csup.items()}
def close_trace(st, maxsteps=600, frozen=None, skip=frozenset()):
    st=st.clone(); frozen=set(frozen or ()); trace=[]
    for step in range(maxsteps):
        nz=[a for a in st.nz() if a not in skip]
        if not nz: return st,True,trace,frozen
        nz.sort(key=lambda a: SUPN[a]); moved=False
        for a in nz:
            best=None
            for X in AVARS[a]:
                if X not in FREE or X in frozen: continue
                for rt in roots(st,a,X):
                    if st.fv.get(X,0)==rt: continue
                    g=st.clone().set_free({X:rt})
                    if g.av[a]!=0: continue
                    k=(len(g.nz()),len(g.fails))
                    if best is None or k<best[0]: best=(k,g,X,rt)
            if best is None: continue
            k,g,X,rt=best
            st=g; frozen.add(X); trace.append((a,X)); moved=True; break
        if not moved: return st,False,trace,frozen
    return st,False,trace,frozen
if __name__=='__main__':
    BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
    st0=St({}); b=BITS['A'][0]
    st=st0.clone().set_free({b:1})
    out,ok,tr,fr=close_trace(st,frozen=set(ALL))
    print('closed=%s score=%d  chain length=%d'%(ok,out.score(),len(tr)))
    print('cascade pin atoms (atom <- free input), in closure order:')
    for a,X in tr: print('   a%-6d <- x_%-6d  (%d equations)'%(a,X,len(ev.eq_terms) and sum(1 for i,(m,sq,tl) in enumerate(ev.eq_terms) if a in [q for c,q in tl])))
    json.dump([[a,X] for a,X in tr],open('chain.json','w'))
