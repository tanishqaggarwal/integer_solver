import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
v=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json'); ad.fwd(v,rounds=6)
def show(a):
    print('--- a%d  (%d eqs)  value=%d'%(a,len(L.atom2eq.get(a,{})),L.evalpoly(L.polys[a],v)))
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(-len(kv[0]),kv[0])):
        vs=' * '.join('x%d[%s]'%(u,('%d'%v[u]) if v[u].bit_length()<25 else str(v[u])[:9]+'..') for u in m)
        print('     %-70s coeff %s'%(vs, c if abs(c)<10**12 else str(c)[:14]+'..(%dd)'%len(str(abs(c)))))
for a in [1618,688,31670,3576,3578,31672,2423,29539,33929,26731,7930,21617,33796,25676,40562,33792,40623]:
    show(a)
