"""Dual formulation: pick an equation set S >= FAILS; every variable whose entire equation
footprint lies inside S is then a FREE knob.  Score = |S| - maxzero(S).  Grow S greedily."""
import sys, pickle, itertools, collections, random, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from dioph2 import solve_int
from opt import atom_linear

d=pickle.load(open('atoms.pkl','rb')); eq_terms=d['eq_terms']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
codes,_=H.load_equations()
FAILS=sorted(H.evaluate(codes,BASE))

# variable -> full equation footprint
varfoot={}
for x in range(NV):
    ats=var_atoms.get(x)
    if not ats: continue
    E=set()
    for a in ats: E |= set(atom2eq.get(a,[]))
    varfoot[x]=E
print('vars with a footprint:',len(varfoot))

def knobs_of(S):
    S=set(S)
    return sorted(x for x,E in varfoot.items() if E and E <= S)

def build_forms(S, V):
    V=set(V)
    forms={}
    for e in sorted(S):
        m,sq,tl=eq_terms[e]
        const=0; coef=collections.defaultdict(int)
        for c,a in tl:
            r=atom_linear(a,V)
            if r is None: return None
            b,cf=r
            const+=c*b
            for v,x in cf.items(): coef[v]+=c*x
        forms[e]=(const,{v:x for v,x in coef.items() if x})
    return forms

def maxzero(S, V, trials=60, seed=0):
    forms=build_forms(S,V)
    if forms is None: return None
    Vl=sorted(V); rnd=random.Random(seed)
    best=(0,(),None)
    Sl=sorted(S)
    for t in range(trials):
        order=Sl[:] if t==0 else rnd.sample(Sl,len(Sl))
        chosen=[]; sol=None
        for e in order:
            tr=chosen+[e]
            M=[[forms[x][1].get(v,0) for v in Vl] for x in tr]
            r=[-forms[x][0] for x in tr]
            s=solve_int(M,r)
            if s is not None: chosen=tr; sol=s
            if len(chosen)>=len(Vl): break
        if len(chosen)>best[0]: best=(len(chosen),tuple(chosen),sol)
    return best[0], best[1], best[2], Vl

def score(S, trials=60):
    V=knobs_of(S)
    if not V: return None
    r=maxzero(S,V,trials=trials)
    if r is None: return None
    k,sub,sol,Vl=r
    return len(S)-k, k, Vl, sub, sol

def realise(Vl,sol):
    v=list(BASE)
    for var,dv in zip(Vl,sol): v[var]=BASE[var]+dv
    return v

if __name__=='__main__':
    t0=time.time()
    S0=set(FAILS)
    r=score(S0); print('S=FAILS:',r[0],'failing; knobs',r[2])
    # candidate equations to add: those sharing atoms with the current S
    cur=set(FAILS)|{9123,18673}
    r=score(cur); print('S=FAILS+{9123,18673}:',r[0],'failing; knobs',r[2],'k=',r[1])
    best=(r[0],frozenset(cur))
    seen={frozenset(cur)}
    beam=[(r[0],frozenset(cur))]
    for depth in range(1,7):
        cands=[]
        for f0,S in beam:
            # equations reachable: those sharing an atom with atoms touched by S
            atoms=set()
            for e in S:
                m,sq,tl=eq_terms[e]
                for c,a in tl: atoms.add(a)
            reach=set()
            for a in atoms: reach |= set(atom2eq.get(a,[]))
            for e in sorted(reach-S):
                ns=frozenset(S|{e})
                if ns in seen: continue
                seen.add(ns)
                rr=score(ns, trials=40)
                if rr is None: continue
                cands.append((rr[0], ns, rr))
        if not cands: break
        cands.sort(key=lambda z:(z[0],len(z[1])))
        beam=[(c[0],c[1]) for c in cands[:5]]
        print(f'depth {depth}: best={cands[0][0]} |S|={len(cands[0][1])} k={cands[0][2][1]} knobs={len(cands[0][2][2])} [{time.time()-t0:.0f}s]',flush=True)
        if cands[0][0]<best[0]:
            best=(cands[0][0],cands[0][1])
            f,k,Vl,sub,sol=cands[0][2]
            v=realise(Vl,sol); ff=H.evaluate(codes,v)
            print('   ACTUAL',len(codes)-len(ff),'/',len(codes),f'({len(ff)} failing)',flush=True)
            if len(ff)<9: H.save_assignment(v,f'pins/s2_{len(ff)}.json')
        pickle.dump(beam,open('pins/search2_beam.pkl','wb'))
