#!/usr/bin/env python3
"""The general compensated-move question, in its sharpest form.

Is there ANY move (however many variables, however long the compensation chain) that
keeps every atom outside the support T exactly zero while changing the 7 cluster
values in a direction the single-variable lattice does not already contain?

Formulation.  Remove T's atoms from the definer map, so the variables they defined
become free.  Propagate derivatives through the remaining DAG (exact forward-mode AD),
which keeps every DEFINER atom outside T identically zero by construction.  What is
left to impose is zero gradient on the CONSTRAINT atoms outside T.  Then:

    new generators exist  <=>  some d with J_constraints . d = 0 has J_T . d != 0
                          <=>  the rows of J_T are NOT in the row space of J_constraints

which is the same frozen/not-frozen test used earlier, now applied to the cluster.
"""
import os, pickle, sys, time
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL
import jdiag as D

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
val = E.load(CL.DEL)
T = CL.all_nonzero(val)
Tset = set(T)
print("T =", T)

# --- definer map with T's atoms removed -------------------------------------
definer = dict(J.definer)
released = [v for v, i in definer.items() if i in Tset]
for v in released:
    del definer[v]
print("variables released by dropping T's definers:", released)

order, cyc = E.topo(definer)
assert not cyc, "cycle after release"
FREE = set(range(E.NV)) - set(definer)
EVI = {}
for v, i in definer.items():
    p = polys[i]
    c = p[(v,)]
    EVI[v] = (c, tuple((k, cc) for k, cc in p.items() if k != (v,)))

CONS = sorted(set(range(len(polys))) - set(definer.values()) - Tset)
print("constraint atoms outside T:", len(CONS))


def dcolumn(z, targets):
    """exact d/dz of each atom in `targets`, over Z, at `val`"""
    der = defaultdict(int)
    der[z] = 1
    for v in order:
        e = EVI.get(v)
        if e is None or v in FREE:
            continue
        c, rest = e
        ds = 0
        for k, cc in rest:
            t = cc
            dt = 0
            for j in k:
                dt = dt * val[j] + t * der[j]
                t *= val[j]
            ds += dt
        if ds:
            der[v] = -ds // c
        else:
            der[v] = 0
    out = {}
    for i in targets:
        s = 0
        for k, c in polys[i].items():
            t = c
            dt = 0
            for j in k:
                dt = dt * val[j] + t * der[j]
                t *= val[j]
            s += dt
        if s:
            out[i] = s
    return out


# knobs: free variables that can reach T at all, plus their compensation partners
seen = set(); q = deque()
for i in T:
    q.extend(E.varsof[i])
lv = set()
while q:
    x = q.popleft()
    if x in seen: continue
    seen.add(x)
    d = definer.get(x)
    if d is None:
        lv.add(x); continue
    for w in E.varsof[d]:
        if w != x and w not in seen: q.append(w)
knobs = sorted(lv)
print("knobs reaching T:", len(knobs))

t0 = time.time()
JT = {}
JC = {}
for jx, z in enumerate(knobs):
    col = dcolumn(z, T + CONS)
    for i, v in col.items():
        if i in Tset:
            JT[(i, jx)] = v
        else:
            JC[(i, jx)] = v
print(f"jacobian built in {time.time()-t0:.0f}s; nnz(J_T)={len(JT)} nnz(J_cons)={len(JC)}")

rowsC = defaultdict(dict)
for (i, jx), v in JC.items():
    rowsC[i][jx] = v
rowsT = defaultdict(dict)
for (i, jx), v in JT.items():
    rowsT[i][jx] = v
print("constraint rows with support:", len(rowsC), " T rows with support:", len(rowsT))

pivC, _ = D.echelon([(r, 0) for r in rowsC.values()])
print(f"rank of the outside-T constraint Jacobian: {len(pivC)} (of {len(knobs)} cols)"
      f"  => kernel dimension {len(knobs)-len(pivC)}")

newdir = 0
for i in T:
    r = rowsT.get(i, {})
    if not r:
        print(f"  a{i}: no gradient in any knob")
        continue
    row, b, dep = D.reduce_row(pivC, r, 0)
    if dep:
        print(f"  a{i}: gradient lies inside the constraint row space => FROZEN "
              f"on the compensated kernel")
    else:
        newdir += 1
        print(f"  a{i}: gradient ESCAPES the constraint row space "
              f"({len(row)} independent cols) => compensated moves DO reach it")
print(f"\ncluster coordinates reachable by compensated moves: {newdir} of {len(T)}")
pickle.dump({'knobs': knobs, 'rowsC': dict(rowsC), 'rowsT': dict(rowsT)},
            open(os.path.join(HERE, 'jpairs2.pkl'), 'wb'))
