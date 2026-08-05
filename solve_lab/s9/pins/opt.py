"""Exact model: pick a set V of variables to move; every affected equation becomes an integer
linear form in the deltas.  Maximise the number of affected equations that can be zeroed."""
import sys, pickle, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from dioph2 import solve_int

d=pickle.load(open('atoms.pkl','rb')); src=d['atom_src']; eq_terms=d['eq_terms']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
roots=pickle.load(open('roots.pkl','rb'))
codes,_=H.load_equations()
FAILS=sorted(H.evaluate(codes,BASE))
freeset=set(x for x in range(NV) if x not in definer)

def atom_linear(a, V):
    """Return (base_value, {v: coef}) if polys[a] is affine in the V-vars, else None."""
    Pp = roots[a] if (a in roots) else polys[a]
    base=0; coef=collections.defaultdict(int)
    for m,c in Pp.items():
        inV=[u for u in m if u in V]
        if len(inV)>1: return None
        t=c
        for u in m:
            if u not in V: t*=BASE[u]
        if inV:
            coef[inV[0]]+=t
            base+= t*BASE[inV[0]]
        else:
            base+=t
    return base, dict(coef)

def model(V):
    """Return (S, forms) where forms[e] = (const, {v:coef}); equation e is zero iff form==0."""
    V=set(V)
    atoms=set()
    for v in V: atoms |= set(var_atoms[v])
    S=set(FAILS)
    for a in atoms: S |= set(atom2eq.get(a,[]))
    forms={}
    for e in sorted(S):
        m,sq,tl=eq_terms[e]
        const=0; coef=collections.defaultdict(int)
        ok=True
        for c,a in tl:
            r=atom_linear(a,V)
            if r is None: ok=False; break
            b,cf=r
            const+=c*b
            for v,x in cf.items(): coef[v]+=c*x
        if not ok: return None
        forms[e]=(const, {v:x for v,x in coef.items() if x})
    return sorted(S), forms

def maxzero(V, verbose=False):
    r=model(V)
    if r is None: return None
    S,forms=r
    Vl=sorted(V)
    best=(0,())
    for k in range(len(Vl),0,-1):
        if k>len(S): continue
        found=None
        for sub in itertools.combinations(S,k):
            M=[[forms[e][1].get(v,0) for v in Vl] for e in sub]
            rr=[-forms[e][0] for e in sub]
            sol=solve_int(M,rr)
            if sol is not None: found=(sub,sol); break
        if found:
            best=(k,found[0],found[1])
            break
    return S, Vl, best

CAND=[9413, 28730, 642, 17325, 1329, 6947, 10903, 23754, 29854, 31864, 33168,
      1844, 21574, 35619, 950, 1613, 9629, 15120, 35531, 10422]

if __name__=='__main__':
    base=[9413,28730]
    r=maxzero(base); S,Vl,best=r
    print(f'V={Vl}  |S|={len(S)}  max zeroed={best[0]}  -> failing={len(S)-best[0]}')
    print('   zeroed:',best[1])
    results=[]
    for extra in range(1,4):
        for combo in itertools.combinations([c for c in CAND if c not in base], extra):
            V=base+list(combo)
            r=maxzero(V)
            if r is None:
                continue
            S,Vl,best=r
            fail=len(S)-best[0]
            results.append((fail,len(S),best[0],tuple(combo)))
        results.sort()
        print(f'\n--- after adding {extra} extra knob(s): best 10 ---')
        for f,ls,k,c in results[:10]:
            print(f'   failing={f}  |S|={ls} zeroed={k}  extra={c}')
    pickle.dump(results, open('pins/opt_results.pkl','wb'))
