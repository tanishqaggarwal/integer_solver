"""Clusters: sets of atoms whose equations OVERLAP, so many can be broken for few equations.

This is how the 39,026 checkpoint really works: it breaks five gate atoms and two checks, but
they all live in the same seven equations, so the price is 7 rather than 5*10.  Generalise it:
enumerate small equation-sets E, take every atom whose equations are contained in E (all of them
can be broken for a total cost of |E|), and test whether the knobs so purchased close the mod-p
system.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHK=set(a for a in range(L.NA) if L.atom_out.get(a) is None)
EQ={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
small=[a for a in range(L.NA) if 0<len(EQ[a])<=7]
print(f"atoms in <=7 equations: {len(small)}")
# group by equation set; then merge sets whose union is still small
byeq=collections.defaultdict(list)
for a in small: byeq[EQ[a]].append(a)
keys=sorted(byeq, key=len)
print(f"distinct equation-sets among them: {len(keys)}")
clusters={}
for i,k1 in enumerate(keys):
    for k2 in keys[i:]:
        u=k1|k2
        if len(u)>7: continue
        atoms=tuple(sorted(a for a in small if EQ[a]<=u))
        if len(u) < clusters.get(atoms,(99,))[0]:
            clusters[atoms]=(len(u),u)
best=sorted(((c,atoms) for atoms,(c,u) in clusters.items()), key=lambda z:(z[0],-len(z[1])))
print("cheapest clusters (cost, #atoms, atoms):")
seen=set()
show=[]
for c,atoms in best:
    if atoms in seen: continue
    seen.add(atoms)
    ng=sum(1 for a in atoms if a not in CHK)
    show.append((c,len(atoms),ng,atoms))
show.sort(key=lambda z:(z[0],-z[2]))
for c,n,ng,atoms in show[:18]:
    print(f"   cost {c} eqs: {n} atoms ({ng} gates) {atoms[:10]}")
json.dump([[c,list(a)] for c,n,ng,a in show[:200]], open(os.path.join(HERE,'data','gmp31.json'),'w'))
