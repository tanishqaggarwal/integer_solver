"""Global combinatorial search: minimize |E(S)| - |S| over atom sets S.
Equivalent dual: maximize f(E) = |A(E)| - |E| where A(E) = {a : eqs(a) subset of E}.
Greedy 1-missing closure growth from many seeds."""
import os, sys, collections, json, time, heapq
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
eq_atoms_idx=[set() for _ in range(L.NEQ)]
for a in range(L.NA):
    for e in EQS[a]: eq_atoms_idx[e].add(a)

def grow(E0, maxsteps=400, allow_zero=8):
    E=set(E0)
    missing={a: len(EQS[a]-E) for a in range(L.NA)}
    # only track atoms with small missing
    best=None
    zeros=0
    for step in range(maxsteps):
        cand=collections.Counter()
        for a,mi in missing.items():
            if mi==1:
                e=next(iter(EQS[a]-E)); cand[e]+=1
        if not cand: break
        e,g=cand.most_common(1)[0]
        gain=g-1
        if gain<0:
            if zeros>=allow_zero: break
            zeros+=1
        E.add(e)
        for a in eq_atoms_idx[e]:
            if a in missing: missing[a]-=1
        A=[a for a in missing if missing[a]==0]
        f=len(A)-len(E)
        if best is None or f>best[0]: best=(f,frozenset(E),tuple(sorted(A)))
    A=[a for a in missing if missing[a]==0]
    f=len(A)-len(E)
    if best is None or f>best[0]: best=(f,frozenset(E),tuple(sorted(A)))
    return best

t0=time.time()
S0=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
E0=set().union(*[EQS[a] for a in S0])
A0=[a for a in range(L.NA) if EQS[a]<=E0]
print('current cluster: |E|=%d atoms fully inside=%s  f=%d'%(len(E0),A0,len(A0)-len(E0)))
b=grow(E0)
print('grown from current: f=%d |E|=%d |A|=%d'%(b[0],len(b[1]),len(b[2])))
print('   A=',b[2][:40])
print('   E=',sorted(b[1]))
# global scan over single-atom seeds, ranked by footprint
order=sorted(range(L.NA),key=lambda a:len(EQS[a]))
res=[]
for i,a in enumerate(order):
    if len(EQS[a])>9: break
    bb=grow(EQS[a],maxsteps=200,allow_zero=6)
    res.append((bb[0],a,len(bb[1]),len(bb[2])))
    if i%200==0: print(i,a,len(EQS[a]),bb[0],f'{time.time()-t0:.0f}s',flush=True)
res.sort(key=lambda t:-t[0])
print('\nTOP seeds by f = |A|-|E|:')
for f,a,ne,na in res[:30]:
    print(f'  seed a{a} fp={len(EQS[a])}  f={f}  |E|={ne} |A|={na}')
json.dump([[f,a,ne,na] for f,a,ne,na in res],open(os.path.join(HERE,'pa_closure.json'),'w'))
print(f'{time.time()-t0:.0f}s')
