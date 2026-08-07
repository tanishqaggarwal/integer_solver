#!/usr/bin/env python3
"""Agent P: structure pass on my own model2.pkl."""
import pickle, sys
from collections import defaultdict, Counter

D = pickle.load(open('/home/user/integer_solver/solve_lab/agentP_work/model2.pkl','rb'))
rows, AP = D['eq_rows'], D['atom_polys']

use = Counter()
for r in rows:
    for c,a in r['row']: use[a]+=1
print("atom slots:", sum(use.values()), "distinct:", len(AP))
print("usage histogram:", sorted(Counter(use.values()).items())[:15])

# classify atoms
def cls(ap):
    vs=set()
    for m in ap: vs.update(m)
    deg=max((len(m) for m in ap), default=0)
    return deg, len(vs)

cl=Counter()
for ap in AP: cl[cls(ap)]+=1
print("atom (deg,nvars):", dict(sorted(cl.items())))

# boolean atoms:  x^2 - x  (any scaling)
def is_bool(ap):
    ks=set(ap)
    if len(ks)!=2: return None
    for m in ks:
        if len(m)==2 and m[0]==m[1]:
            v=m[0]
            if (v,) in ap and ap[m]==-ap[(v,)]: return v
    return None

nb=0; boolvars=set()
for i,ap in enumerate(AP):
    v=is_bool(ap)
    if v is not None: nb+=1; boolvars.add(v)
print("boolean-form atoms:", nb, "distinct bool vars:", len(boolvars))

# pin atoms: x - c
pins=defaultdict(set)
for i,ap in enumerate(AP):
    if len(ap)<=2:
        vs=[m for m in ap if len(m)==1]
        if len(vs)==1 and set(ap)-{vs[0]} <= {()}:
            if ap[vs[0]] in (1,-1):
                c = -ap.get((),0)//ap[vs[0]]
                pins[vs[0][0]].add(c)
print("pin atoms (x-c):", sum(len(v) for v in pins.values()), "vars pinned:", len(pins))
multi=[ (k,v) for k,v in pins.items() if len(v)>1 ]
print("vars with >1 distinct pin value:", len(multi), multi[:5])
allv=set()
for ap in AP:
    for m in ap: allv.update(m)
print("distinct variables appearing:", len(allv), "max idx", max(allv))
