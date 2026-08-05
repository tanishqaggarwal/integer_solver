"""Beam search over knob-variable sets, maximising (#zeroable equations - #affected equations)."""
import sys, pickle, itertools, collections, random, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from dioph2 import solve_int
from opt import atom_linear, model

d=pickle.load(open('atoms.pkl','rb')); eq_terms=d['eq_terms']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
codes,_=H.load_equations()
FAILS=sorted(H.evaluate(codes,BASE))

def greedy_maxzero(S, forms, Vl, trials=40, seed=0):
    rnd=random.Random(seed)
    best=(0,(),None)
    order0=list(S)
    for t in range(trials):
        order=order0[:] if t==0 else rnd.sample(order0,len(order0))
        chosen=[]; sol=None
        for e in order:
            trial=chosen+[e]
            M=[[forms[x][1].get(v,0) for v in Vl] for x in trial]
            r=[-forms[x][0] for x in trial]
            s=solve_int(M,r)
            if s is not None:
                chosen=trial; sol=s
            if len(chosen)>=len(Vl): break
        if len(chosen)>best[0]: best=(len(chosen),tuple(chosen),sol)
    return best

def evaluate_V(V, trials=40):
    r=model(V)
    if r is None: return None
    S,forms=r
    Vl=sorted(V)
    k,sub,sol=greedy_maxzero(S,forms,Vl,trials=trials)
    return len(S)-k, S, Vl, sub, sol

def pool_from(V):
    """Variables occurring in atoms of the equations currently affected."""
    r=model(V)
    if r is None: return []
    S,_=r
    atoms=set()
    for e in S:
        m,sq,tl=eq_terms[e]
        for c,a in tl: atoms.add(a)
    vs=set()
    for a in atoms: vs |= set(u for m_ in polys[a] for u in m_)
    return sorted(vs)

if __name__=='__main__':
    t0=time.time()
    beam=[(9, tuple(sorted([9413,28730,642,29854,31864])))]
    seen=set(b[1] for b in beam)
    BW=6
    for depth in range(1,9):
        cands=[]
        for score,V in beam:
            pool=pool_from(V)
            for x in pool:
                if x in V: continue
                nv=tuple(sorted(set(V)|{x}))
                if nv in seen: continue
                seen.add(nv)
                r=evaluate_V(nv, trials=25)
                if r is None: continue
                f=r[0]
                cands.append((f,nv))
        cands.sort()
        if not cands: break
        beam=cands[:BW]
        print(f'depth {depth}: best failing={beam[0][0]}  V={beam[0][1]}  [{time.time()-t0:.0f}s]',flush=True)
        for f,V in beam[:BW]: print(f'    {f}  {V}')
        pickle.dump(beam, open('pins/beam.pkl','wb'))
