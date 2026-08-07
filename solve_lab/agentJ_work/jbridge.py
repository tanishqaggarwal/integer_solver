#!/usr/bin/env python3
"""The a23328 bridge.

a23328 = ((x_4432)-(x_19964))-(x_28730) is the ONLY atom that is not equation-disjoint
from the deliverable's cluster: 10 of its 11 equations lie inside the cluster's 12, and
it brings exactly 1 new row (eq8680).  It is not a confined knob, so breaking it (=
releasing x_4432) disturbs constraints outside the support.  Question: are those
disturbances repairable, and what is the exact optimum of the resulting support?
"""
import os, sys, pickle
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL
import jdist as DI

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
val0 = E.load(CL.DEL)
DELT = CL.all_nonzero(val0)
print("deliverable support:", DELT)
print("a23328 equations:", sorted(DI.A2E[23328]))
Rdel = set()
for j in DELT:
    Rdel |= DI.A2E[j]
print("cluster rows:", sorted(Rdel))
print("new rows a23328 brings:", sorted(DI.A2E[23328] - Rdel))

# --- what does releasing x_4432 disturb, exactly? ---------------------------
definer = dict(J.definer)
print("\na23328 is the definer of:", [v for v, i in definer.items() if i == 23328])
rel = [v for v, i in definer.items() if i in set(DELT) | {23328}]
for v in rel:
    del definer[v]
order, cyc = E.topo(definer)
assert not cyc
FREE = set(range(E.NV)) - set(definer)
EVI = {}
for v, i in definer.items():
    p = polys[i]
    EVI[v] = (p[(v,)], tuple((k, cc) for k, cc in p.items() if k != (v,)))
uses = defaultdict(list)
for v, i in definer.items():
    for w in E.varsof[i]:
        if w != v:
            uses[w].append(v)
pos = {v: k for k, v in enumerate(order)}


def fwd(val, seed=None):
    if seed is None:
        it = order
    else:
        dirty = set(); q = deque([seed])
        while q:
            x = q.popleft()
            for w in uses[x]:
                if w not in dirty and w not in FREE:
                    dirty.add(w); q.append(w)
        it = sorted(dirty, key=lambda z: pos[z])
    for v in it:
        e = EVI.get(v)
        if e is None or v in FREE:
            continue
        c, rest = e
        s = 0
        for k, cc in rest:
            t = cc
            for j in k:
                t *= val[j]
            s += t
        val[v] = -s // c
    return val


assert fwd(list(val0)) == val0, "deliverable not a fixed point of the released DAG"
print("released variables:", sorted(rel))

for d in (1, 2, 1000003):
    v1 = list(val0)
    v1[4432] += d
    fwd(v1, seed=4432)
    nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
    extra = sorted(set(nz) - set(DELT) - {23328})
    print(f"  delta={d}: nonzero atoms {len(nz)}, outside cluster+a23328: {extra}")

# --- the resulting support and its exact bounds -----------------------------
v1 = list(val0); v1[4432] += 1
fwd(v1, seed=4432)
NZ = sorted(i for i in range(len(polys)) if CL.atomval(i, v1) != 0)
print("\nsupport after breaking the bridge:", NZ)
a, r = DI.stats(NZ)[0], DI.stats(NZ)[1]
print(f"  alone={a} |R|={r}  => weight >= {max(a, 0)}   (deliverable weight 7)")
for T in ([*DELT, 23328], NZ):
    a, r, t = DI.stats(T)
    print(f"  support {sorted(set(T))[:4]}... |T|={t} |R|={r} alone={a} "
          f"generic |R|-|T|+1={r-t+1}")
pickle.dump({'NZ': NZ}, open(os.path.join(HERE, 'jbridge.pkl'), 'wb'))
