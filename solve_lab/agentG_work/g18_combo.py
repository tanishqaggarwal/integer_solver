import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gred, gGclose
from gsym import *
SRC='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v0=L.load(SRC); ad.fwd(v0,rounds=6)
combos=[[],[2081],[24601],[2081,24601]]
if len(sys.argv)>1: combos=[[int(x) for x in c.split(',') if x] for c in sys.argv[1:]]
for FL in combos:
    v=list(v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    sc=L.NEQ-len(L.failing_eqs(L.all_atom_values(v)))
    S=gGclose.closure(v)
    r=gred.reduce_state(v,S)
    res=[(a,(g%P if isinstance(g,int) else 'POLY%d'%len(g))) for a,g in r['res'] if (isinstance(g,int) and g%P) or not isinstance(g,int)]
    print('flip %-14s raw=%d |S|=%d lin=%d rank=%d ninc=%d nfree=%d nonlin=%d nzc=%d residual=%s'
          %(FL,sc,len(S),r['nlin'],r['rank'],r['ninc'],r['nfree'],r['nnon'],len(r['nzc']),
            [(a,str(g)[:10]) for a,g in res]), flush=True)
