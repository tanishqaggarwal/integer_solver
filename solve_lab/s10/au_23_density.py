import os, sys, collections, itertools, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E12=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
eqs={a:frozenset(L.atom2eq[a]) for a in range(L.NA)}
# check identical equation-sets
byeq=collections.defaultdict(list)
for a in range(L.NA): byeq[eqs[a]].append(a)
dupes={k:v for k,v in byeq.items() if len(v)>1}
print('groups of atoms with IDENTICAL equation sets:', len(dupes), ' total atoms in them:', sum(len(v) for v in dupes.values()))
big=sorted(dupes.values(), key=lambda v:-len(v))[:6]
for g in big: print('   size',len(g),'atoms',g[:8],'neq',len(eqs[g[0]]))

# --- minimise |E(N)| - |N| over N containing the seven, greedy + local search ---
base=frozenset().union(*[eqs[a] for a in SEVEN])
print('\nE(seven) =',len(base), sorted(base))
N=set(SEVEN); E=set(base)
POOL=[a for a in range(L.NA) if a not in N and len(eqs[a]-E)<=3 and len(eqs[a])<=40]
print('pool of atoms adding <=3 new equations:', len(POOL))
improved=True
while improved:
    improved=False
    cands=sorted(((len(eqs[a]-E),a) for a in range(L.NA) if a not in N), key=lambda t:t[0])
    for d,a in cands:
        if d<=0:
            N.add(a); E|=eqs[a]; improved=True
    if improved: continue
    # allow d==1 only if it's neutral: no gain; stop
    break
print(f'after adding all zero-new-equation atoms: n={len(N)} |E|={len(E)}  |E|-n={len(E)-len(N)}')
print('   added:', sorted(N-set(SEVEN)))
# now full greedy allowing d==1 (neutral) and d==0
N2=set(SEVEN); E2=set(base)
order=[]
while True:
    best=None
    for a in range(L.NA):
        if a in N2: continue
        d=len(eqs[a]-E2)
        if best is None or d<best[0]: best=(d,a)
    if best is None or best[0]>1: break
    d,a=best; N2.add(a); E2|=eqs[a]; order.append((a,d))
    if len(N2)>400: break
print(f'greedy with d<=1: n={len(N2)} |E|={len(E2)} |E|-n={len(E2)-len(N2)}  (added {len(order)})')
print('   first 25 added:', order[:25])
