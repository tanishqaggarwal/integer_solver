import sys, os, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
FREE=[u for u in range(L.NVARS) if u not in L.definer]
T=[11150,25739,37758]
ALLCHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
CH=set(ALLCHK)
mov=[]
for u in FREE:
    v2=list(v); L.ripple(v2,{u:v[u]+1})
    if any((v2[t]-v[t])%P for t in T): mov.append(u)
print("movers:",len(mov))
kinds=collections.Counter('bool' if v[u] in (0,1) else 'other' for u in mov)
print("  by current value:",dict(kinds))
# affine test + collateral count for each mover
rows=[]
for u in mov:
    d=[]
    ok=True
    for s in (1,2):
        v2=list(v); L.ripple(v2,{u:v[u]+s})
        d.append([ (v2[t]-v[t]) for t in T ])
        if s==1:
            diff=[w for w in range(L.NVARS) if v2[w]!=v[w]]
            coll=set()
            for w in diff: coll |= {a for a in L.var_atoms[w] if a in CH}
            nbroke=sum(1 for a in coll if L.evalpoly(L.polys[a],v2)!=0)
            base=sum(1 for a in coll if L.evalpoly(L.polys[a],v)!=0)
    aff = all(d[1][i]==2*d[0][i] for i in range(3))
    rows.append((u, aff, v[u], [x%P for x in d[0]], nbroke-base))
aff=[r for r in rows if r[1]]
print("  affine movers:",len(aff))
aff.sort(key=lambda r:r[4])
print("  cheapest affine movers (u, value, extra broken checks):")
for r in aff[:15]: print(f"     x{r[0]} val={str(r[2])[:14]} extra_broken={r[4]}")
import json
json.dump([[r[0],r[1],r[4]]+[str(x) for x in r[3]] for r in rows], open(os.path.join(HERE,'data','modp2.json'),'w'))
