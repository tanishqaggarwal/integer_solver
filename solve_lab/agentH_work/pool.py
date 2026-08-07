"""Candidate detach pool: gate outputs feeding the residual region."""
import ev, json, pickle
from collections import defaultdict
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
REGION_ATOMS=[22229,22230,22231,35758,35759,35760,35761,35762,37887]
FREE0=set(ev.F['free0'])
pool=set()
for a in REGION_ATOMS:
    for v in ev.atom_vars[a]:
        if v not in FREE0: pool.add(v)          # defined => detachable gate output
# one level up: the definers' inputs that are themselves defined
lvl2=set()
for v in list(pool):
    d=ev.F['definer'][v]
    if d>=0:
        for u in ev.atom_vars[d]:
            if u!=v and u not in FREE0: lvl2.add(u)
pool|=lvl2
WIT=[642,28730,29854,31864]
for v in WIT:
    assert v not in FREE0, v
    pool.add(v)
pool=sorted(pool)
print('detach pool size:',len(pool))
print('pool:',pool)
print('witness set inside pool:',all(v in pool for v in WIT))
json.dump(pool,open('pool.json','w'))
