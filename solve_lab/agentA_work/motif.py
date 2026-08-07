"""Structural motif census: is there a ~256-fold repeated sub-circuit (a doubling ladder)?
Shape of an atom = sorted multiset of (|coeff| class, monomial degree)."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
def cls(c):
    a=abs(c)
    if a<1000: return 's%d'%a
    if a==P: return 'p'
    return 'big%d'%(len(str(a))//10)
shapes=collections.Counter()
for a in range(L.NA):
    sh=tuple(sorted((cls(c),len(m)) for m,c in L.polys[a].items()))
    shapes[sh]+=1
print('distinct atom shapes: %d over %d atoms'%(len(shapes),L.NA))
print('top shape multiplicities:')
for sh,c in shapes.most_common(15):
    print('   %-6d  %s'%(c,str(sh)[:110]))
# how many shape-classes have multiplicity in [200,320]?
band=[(sh,c) for sh,c in shapes.items() if 200<=c<=320]
print('\nshape classes with multiplicity in [200,320] (ladder-sized): %d'%len(band))
for sh,c in sorted(band,key=lambda t:-t[1])[:10]:
    print('   %-6d %s'%(c,str(sh)[:110]))
# equation-shape census
esh=collections.Counter()
for i in range(L.NEQ):
    m,sq,co=L.eq_atoms[i]
    esh[(len(co),sq)]+=1
print('\nequation (n_atoms, is_square) histogram (top 12):',esh.most_common(12))
# degree-4 atoms: how many, and their shapes
d4=[a for a in range(L.NA) if any(len(m)==4 for m in L.polys[a])]
print('\ndegree-4 atoms: %d'%len(d4))
