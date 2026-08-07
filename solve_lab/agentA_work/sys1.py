import sys, json, collections; sys.path.insert(0,'.')
import env, lib as L
from math import gcd
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
S=env.SEVEN
import ahandles as HH
solo,gran=HH.build(v)
H=[a for a in solo if gran[a] and a not in S]
print('handle atoms usable: %d  (gran1=%d, granp=%d)'%(len(H),
      sum(1 for a in H if gran[a]==1), sum(1 for a in H if gran[a]==P)))
R=set()
for a in S+H: R |= set(L.atom2eq[a])
print('equations touched by S u H: %d of %d'%(len(R), L.NEQ))
E=sorted(set(e for a in S for e in L.atom2eq[a]))
print('E (12):',E)
# for each eq in E: which atoms, and which have handles
for e in E:
    m,sq,co=L.eq_atoms[e]
    parts=[]
    for a,c in sorted(co.items()):
        tag='S' if a in S else ('H1' if gran.get(a)==1 else ('Hp' if gran.get(a)==P else '.'))
        parts.append('%d:%+d%s'%(a,c,tag))
    print(' eq%-6d %s'%(e,' '.join(parts)))
