import os, sys, json, pickle, time, collections, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
SYMS = json.load(open('closed_nonbool.json'))
v = L.load(src); ad.fwd(v, rounds=6)
print('base score', L.NEQ-len(L.failing_eqs(L.all_atom_values(v))), 'syms', len(SYMS))
t0=time.time()
val = gsym.build(v, SYMS, cap=None)
n=len(SYMS)
rows=[];nzc=[]
for a in gsym.check_atoms():
    f=gsym.evalpoly_sym(a,val,n,None)
    if isinstance(f,int):
        if f%P: nzc.append((a,f%P))
    else: rows.append((a,f))
print('symbolic checks %d, nonzero-const %d, %.1fs'%(len(rows),len(nzc),time.time()-t0))
print('deg hist', dict(sorted(collections.Counter(gsym.deg(f) for a,f in rows).items())))
print('total terms', sum(len(f) for a,f in rows))
print('nonzero constants (unreachable):', nzc[:20])
pickle.dump({'syms':SYMS,'rows':rows,'nzc':nzc}, open('sys112.pkl','wb'))
# validation at random point
random.seed(11); n=len(SYMS)
def ev(f,pt):
    if isinstance(f,int): return f%P
    s=0
    for m,c in f.items():
        t=c
        for i,e in enumerate(m):
            if e: t=t*pow(pt[i],e,P)%P
        s=(s+t)%P
    return s
for trial in range(2):
    pt=[random.randrange(P) for _ in range(n)]
    w=list(v)
    for i,u in enumerate(SYMS): w[u]=pt[i]
    ad.fwd(w,rounds=10)
    av=L.all_atom_values(w); D=dict(rows); bad=0
    for a in gsym.check_atoms():
        pred=ev(D[a],pt) if a in D else (dict(nzc).get(a,0))
        if pred!=av[a]%P:
            bad+=1
            if bad<6: print('  MISMATCH a%d'%a)
    print('validation trial %d mismatches=%d'%(trial,bad))
