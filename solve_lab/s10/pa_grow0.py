"""Exact level-k compensator analysis around E0 (the delivered cluster's 12 equations).
Question: is there an equation set X such that |A(E0 u X)| - |E0 u X| > |A(E0)| - |E0| = -4 ?
A(E) counts PRIMITIVE atoms only (bundles cost +1 to c so they never help)."""
import os, sys, collections, json, itertools, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
PRIM=[a for a in range(L.NA) if len(EQS[a])>1]
S0=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
E0=frozenset().union(*[EQS[a] for a in S0])
miss={b:(EQS[b]-E0) for b in PRIM}
h=collections.Counter(len(m) for m in miss.values())
print('E0',sorted(E0))
print('missing hist:',sorted(h.items())[:8])
for k in range(0,5):
    lst=[b for b,m in miss.items() if len(m)==k]
    print(f'  missing={k}: n={len(lst)} {lst[:25]}')
# exact max-closure over E >= E0 restricted to atoms with missing<=6, by min-cut style
# small enough: use greedy + exhaustive over extra-equation sets of size<=3 from the
# equations that appear in some atom's missing set of size <=3
pool=set()
for b,m in miss.items():
    if 1<=len(m)<=3: pool|=m
pool=sorted(pool)
print('candidate extra equations (from atoms missing<=3):',len(pool))
small={b:m for b,m in miss.items() if 1<=len(m)<=3}
def f_of(X):
    X=set(X)
    n=sum(1 for b,m in small.items() if m<=X)
    return len(S0)+1+n-len(E0)-len(X)   # |A(E0)|=8 (S0 + a22231)
print('f(E0) =',f_of(()))
best=[(f_of(()),())]
t0=time.time()
for k in (1,2):
    bb=None
    for X in itertools.combinations(pool,k):
        g=f_of(X)
        if bb is None or g>bb[0]: bb=(g,X)
    print(f'  best |X|={k}: f={bb[0]} X={bb[1]}  ({time.time()-t0:.0f}s)',flush=True)
    best.append(bb)
# size 3 restricted to the 60 most useful equations
cnt=collections.Counter()
for b,m in small.items():
    for e in m: cnt[e]+=1
top=[e for e,_ in cnt.most_common(70)]
bb=None
for X in itertools.combinations(top,3):
    g=f_of(X)
    if bb is None or g>bb[0]: bb=(g,X)
print(f'  best |X|=3 (top70): f={bb[0]} X={bb[1]}  ({time.time()-t0:.0f}s)')
# also: greedy unlimited growth
E=set(E0); n=8
for step in range(40):
    cand=collections.Counter()
    for b,m in miss.items():
        r=m-E
        if len(r)==1: cand[next(iter(r))]+=1
    if not cand: break
    e,g=cand.most_common(1)[0]
    E.add(e)
    n=sum(1 for b,m in miss.items() if m<=E)+len(S0)+1-len([b for b in S0 if False])
    n=sum(1 for b in PRIM if EQS[b]<=E)
    print(f'   greedy step {step}: +eq{e} gains {g} atoms -> |E|={len(E)} |A|={n} f={n-len(E)}')
    if n-len(E)>-4: print('    *** BEATS -4')
