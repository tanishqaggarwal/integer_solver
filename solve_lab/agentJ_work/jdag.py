#!/usr/bin/env python3
"""Build the circuit: orient atoms as definitions, find free inputs, topo order."""
import os, pickle, sys, time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
polys, defcands = P['polys'], P['defcands']
NA = len(polys)
NV = 38748

# --- classify atoms by exact monomial signature -------------------------------
def sig(p):
    d = Counter(len(k) for k in p)
    return (d[0], d[1], d[2])

cls = Counter(sig(p) for p in polys)
print("signature (const,lin,quad) histogram:", cls.most_common(20))

# atoms with no def candidate
nod = [i for i in range(NA) if not defcands[i]]
print("no-def atoms:", len(nod))
sub = Counter(sig(polys[i]) for i in nod)
print("  their signatures:", sub.most_common(20))

# variable occurrence
occ = defaultdict(list)
for i, p in enumerate(polys):
    vs = set()
    for k in p:
        vs.update(k)
    for v in vs:
        occ[v].append(i)
print("vars occurring:", len(occ), "of", NV)
print("occ histogram:", Counter(len(v) for v in occ.values()).most_common(12))

# how many vars are the *unique* def candidate of some atom
defby = defaultdict(list)
for i in range(NA):
    for v in defcands[i]:
        defby[v].append(i)
print("vars that are a def candidate somewhere:", len(defby))

pickle.dump({'occ': dict(occ), 'defby': dict(defby)},
            open(os.path.join(HERE, 'jocc.pkl'), 'wb'))
