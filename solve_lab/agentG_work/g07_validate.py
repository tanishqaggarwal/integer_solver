import os, sys, json, pickle, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('sys17.pkl','rb')); SYMS=d['syms']; rows=d['rows']
src='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v0=L.load(src); ad.fwd(v0,rounds=6)
n=len(SYMS)
random.seed(7)
def ev(f, pt):
    if isinstance(f,int): return f%P
    s=0
    for m,c in f.items():
        t=c
        for i,e in enumerate(m):
            if e: t=t*pow(pt[i],e,P)%P
        s=(s+t)%P
    return s
for trial in range(3):
    pt=[random.randrange(P) for _ in range(n)] if trial else [v0[u]%P for u in SYMS]
    v=list(v0)
    for i,u in enumerate(SYMS): v[u]=pt[i]
    ad.fwd(v,rounds=8)
    av=L.all_atom_values(v)
    bad=0; checked=0
    symset=set(a for a,f in rows)
    for a in gsym.check_atoms():
        pred = ev(dict(rows)[a],pt) if a in symset else 0
        act = av[a]%P
        checked+=1
        if pred!=act:
            bad+=1
            if bad<6: print('  MISMATCH a%d pred=%d act=%d'%(a,pred,act))
    print('trial %d: checks=%d mismatches=%d' % (trial,checked,bad))
