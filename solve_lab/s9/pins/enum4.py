import sys, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from opt import model
from dioph2 import solve_int
import harness as H
codes,_=H.load_equations()
V=[642,9413,28730,29854,31864]
S,forms=model(V); Vl=sorted(V)
best=None; ok=0
for sub in itertools.combinations(S,4):
    M=[[forms[e][1].get(v,0) for v in Vl] for e in sub]
    r=[-forms[e][0] for e in sub]
    s=solve_int(M,r)
    if s is None: continue
    ok+=1
    v=list(BASE)
    for var,dv in zip(Vl,s): v[var]=BASE[var]+dv
    f=H.evaluate(codes,v)
    if best is None or len(f)<best[0]:
        best=(len(f),sub,list(v),f); print('new best',len(f),sub,f,flush=True)
print('integer-solvable 4-subsets:',ok,'of',len(list(itertools.combinations(S,4))))
print('BEST',best[0],best[1],best[3])
if best[0]<9: H.save_assignment(best[2],'pins/enum_%d.json'%best[0])
