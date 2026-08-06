"""Mod-p repair iteration.

Mod p there is no divisibility obstruction at all: every "p | X" condition becomes the plain
equation "X == 0 in GF(p)", and solving a check for a free input never fails on divisibility --
it fails only if the coefficient vanishes mod p.  That is exactly what blocked every integer
repair, so the iteration should behave completely differently here.

Loop: forward-evaluate mod p, take the failing checks, solve each for a free input occurring
linearly in it, apply, repeat.
"""
import sys, os, json, time, collections, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp, solvep, coefp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)

def failing(v): return [a for a in CHK if evalp(L.polys[a],v)]

def repair(v, rounds=200, batch=True, seed=0, verbose=True):
    rnd=random.Random(seed)
    best=None
    hist=[]
    for it in range(rounds):
        forwardp(v)
        F=failing(v)
        if best is None or len(F)<best[0]: best=(len(F),list(v))
        if verbose and (it<8 or it%10==0):
            print(f"  it{it}: failing checks mod p = {len(F)}  {F[:12]}", flush=True)
        if not F: return v, []
        hist.append(len(F))
        if len(hist)>25 and len(set(hist[-25:]))==1: 
            if verbose: print("  stalled"); break
        seeds={}
        order=list(F); rnd.shuffle(order)
        for a in order:
            cands=[t for t in L.avars[a] if t in FREE and t not in seeds]
            rnd.shuffle(cands)
            for t in cands:
                x=solvep(a,t,v)
                if x is not None and x!=v[t]:
                    seeds[t]=x; break
            if not batch and seeds: break
        if not seeds:
            if verbose: print("  no free input can move any failing check"); break
        for t,x in seeds.items(): v[t]=x
    forwardp(v)
    return v, failing(v)

if __name__=='__main__':
    base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
    for seed in range(int(sys.argv[1]) if len(sys.argv)>1 else 1):
        v=list(base)
        t0=time.time()
        v,F=repair(v, seed=seed, verbose=(seed==0))
        print(f"seed {seed}: final failing checks mod p = {len(F)}  {F[:14]}  ({time.time()-t0:.0f}s)",
              flush=True)
        if not F:
            json.dump([int(x) for x in v], open(os.path.join(HERE,'data','gmp3_solved.json'),'w'))
            print("  *** ALL CHECKS ZERO MOD P -- saved data/gmp3_solved.json")
            break
