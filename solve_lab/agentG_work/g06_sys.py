import os, sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
SYMS = json.load(open('supp8.json'))['nonbool'] if len(sys.argv)<3 else [int(x) for x in sys.argv[2].split(',')]
CAP = int(sys.argv[3]) if len(sys.argv)>3 else 8
v = L.load(src); ad.fwd(v, rounds=6)
print('base score', L.NEQ-len(L.failing_eqs(L.all_atom_values(v))))
print('symbols(%d): %s' % (len(SYMS),SYMS))
t0=time.time()
val = gsym.build(v, SYMS, cap=CAP)
n=len(SYMS)
rows=[]; nzc=[]
for a in gsym.check_atoms():
    f = gsym.evalpoly_sym(a, val, n, CAP)
    if isinstance(f,int):
        if f%P: nzc.append(a)
        continue
    rows.append((a,f))
print('symbolic checks: %d ; nonzero-constant checks: %d ; %.1fs' % (len(rows),len(nzc),time.time()-t0))
import collections
dd=collections.Counter(gsym.deg(f) for a,f in rows)
print('degree histogram:', dict(sorted(dd.items())))
print('total terms', sum(len(f) for a,f in rows))
# which symbolic checks are currently ZERO (constraints to preserve) vs NONZERO (to fix)
av=L.all_atom_values(v)
cur_nz=[a for a,f in rows if av[a]%P]
print('symbolic checks currently nonzero mod p:', cur_nz)
print('nonzero-const checks (unfixable by these syms):', nzc)
pickle.dump({'syms':SYMS,'rows':rows,'nzc':nzc,'cap':CAP}, open('sys17.pkl','wb'))
