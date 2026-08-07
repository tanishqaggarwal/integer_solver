"""Level-k compensator search around the delivered cluster's 12 equations.
For each atom b compute missing(b) = |eqs(b) \ E0|.  Then look for small equation
sets X with |{b : eqs(b) subset of E0 u X}| - |X| > current, i.e. net gain."""
import os, sys, collections, json, itertools, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
S0=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
E0=frozenset().union(*[EQS[a] for a in S0])
print('E0 =',sorted(E0),len(E0))
miss={}
for b in range(L.NA):
    m=EQS[b]-E0
    if len(m)<=4: miss[b]=m
h=collections.Counter(len(m) for m in miss.values())
print('missing-count histogram (atoms with <=4 missing):',sorted(h.items()))
for k in range(5):
    print(f'  missing={k}: {[b for b,m in miss.items() if len(m)==k]}')
# atoms with missing 0 -> free compensators
base=[b for b,m in miss.items() if len(m)==0]
print('base A(E0) =',base,'  f =',len(base)-len(E0))
# level-1..3: pick X of size k from the union of small missing sets
cands=collections.Counter()
for b,m in miss.items():
    if 1<=len(m)<=3:
        for e in m: cands[e]+=1
print('candidate extra equations:',len(cands))
best=[]
pool=sorted(cands)
def gain(X):
    E=E0|set(X)
    A=[b for b,m in miss.items() if m<=set(X)]
    return len(A)-len(E), A
g0,_=gain(())
print('f(E0) =',g0)
for k in (1,2,3):
    bb=None
    t0=time.time()
    for X in itertools.combinations(pool,k):
        g,A=gain(X)
        if bb is None or g>bb[0]: bb=(g,X,A)
    print(f'best with {k} extra equations: f={bb[0]} X={bb[1]} |A|={len(bb[2])} A={bb[2]}  ({time.time()-t0:.0f}s)',flush=True)
    if k==2 and len(pool)>150: break
