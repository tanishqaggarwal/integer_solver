"""Hill-climb over a curated knob pool (add/remove/swap one variable)."""
import sys, pickle, itertools, random, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from search import evaluate_V
import harness as H
codes,_=H.load_equations()

POOL=[9413,28730,642,1329,6947,17325,10903,23754,29854,31864,33168,1844,21574,
      35619,950,1613,9629,6090,15120,29305,2892,18253,35531,10422,28355,37720,
      23642,11099,22526,37413,23822,34868,7945,17065,7075,14199,23318,12143,
      27500,30108,13502,8731,9118,8976,3629,16495,37254,15324,6247,1956,34600]

def score(V):
    r=evaluate_V(tuple(sorted(V)), trials=30)
    return (10**9,None) if r is None else (r[0], r)

if __name__=='__main__':
    t0=time.time()
    cur=set([9413,28730,642,29854,31864])
    curf,curr=score(cur)
    print('start',sorted(cur),'failing',curf,flush=True)
    memo={tuple(sorted(cur)):curf}
    improved=True; it=0
    while improved and it<12:
        improved=False; it+=1
        moves=[]
        for x in POOL:
            if x in cur: moves.append(('rm',x))
            else: moves.append(('add',x))
        random.Random(it).shuffle(moves)
        for op,x in moves:
            nv=set(cur); nv.discard(x) if op=='rm' else nv.add(x)
            if not nv: continue
            key=tuple(sorted(nv))
            if key in memo: continue
            f,r=score(nv); memo[key]=f
            if f<curf:
                cur, curf, curr = nv, f, r
                print(f'  it{it} {op} x_{x} -> failing {f}  V={sorted(cur)}  [{time.time()-t0:.0f}s]',flush=True)
                improved=True
                # realise + verify
                S,Vl,sub,sol=r[1],r[2],r[3],r[4]
                v=list(BASE)
                for var,dv in zip(Vl,sol): v[var]=BASE[var]+dv
                ff=H.evaluate(codes,v)
                print(f'     ACTUAL {len(codes)-len(ff)}/{len(codes)} ({len(ff)} failing)',flush=True)
                if len(ff)<=9: H.save_assignment(v,f'pins/hill_{len(ff)}.json')
                break
    print('final',sorted(cur),curf)
