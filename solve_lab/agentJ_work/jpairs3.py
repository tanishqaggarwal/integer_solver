#!/usr/bin/env python3
"""Extract the compensated generators explicitly, verify them over Z, union them with
the single-variable lattice, and re-run the exhaustive HNF six-subset search."""
import os, pickle, sys, time, itertools
from fractions import Fraction
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL
import jpairs as P1

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
val0 = E.load(CL.DEL)
T = CL.all_nonzero(val0)
Tset = set(T)
a0 = [CL.atomval(i, val0) for i in T]

definer = dict(J.definer)
for v in [v for v, i in definer.items() if i in Tset]:
    del definer[v]
order, cyc = E.topo(definer)
assert not cyc
FREE = set(range(E.NV)) - set(definer)
EVI = {}
for v, i in definer.items():
    p = polys[i]
    EVI[v] = (p[(v,)], tuple((k, cc) for k, cc in p.items() if k != (v,)))


def forward(val):
    for v in order:
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


# sanity: x* must be a fixed point of the released DAG
chk = forward(list(val0))
assert chk == val0, "deliverable is not a fixed point of the released DAG"

D2 = pickle.load(open(os.path.join(HERE, 'jpairs2.pkl'), 'rb'))
knobs, rowsC = D2['knobs'], D2['rowsC']
n = len(knobs)
print("knobs:", [f"x_{z}" for z in knobs])

basis = P1.int_nullspace(list(rowsC.values()), n)
print("integer kernel dimension:", len(basis))

gens, gnames, seen = [], [], set()
for bi, d in enumerate(basis):
    for scale in (1, 2):
        v1 = list(val0)
        for k, c in enumerate(d):
            if c:
                v1[knobs[k]] += c * scale
        forward(v1)
        nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
        if set(nz) <= Tset:
            dv = tuple(CL.atomval(i, v1) - a0[k] for k, i in enumerate(T))
            if any(dv) and dv not in seen:
                seen.add(dv)
                gens.append(list(dv))
                sup = [knobs[k] for k, c in enumerate(d) if c]
                gnames.append(f"ker{bi}x{scale}[{len(sup)}v:{sup[:4]}]")
        else:
            if scale == 1:
                print(f"  ker{bi}: leaks to {len(set(nz)-Tset)} atoms outside T (second order)")
        break

print(f"\nVERIFIED compensated generators: {len(gens)}")
for n_, g in zip(gnames, gens):
    print("   ", n_, [x if abs(x) < 10**10 else f"~{len(str(abs(x)))}d" for x in g])

# union with the previously found single-variable generators
prev = pickle.load(open(os.path.join(HERE, 'jgen3.pkl'), 'rb'))
allg = [list(g) for g in prev['gens']]
alln = list(prev['gnames'])
for g, nm in zip(gens, gnames):
    if tuple(g) not in {tuple(x) for x in allg}:
        allg.append(g); alln.append(nm)
print(f"\nUNION lattice: {len(allg)} generators")

rows, reqs = prev['rows'], prev['reqs']
nT = len(T)
print("rows:", len(rows))
best = None
for k in range(len(rows), 4, -1):
    found = None; cnt = 0
    for S in itertools.combinations(range(len(rows)), k):
        cnt += 1
        B = [[sum(rows[r][c] * allg[g][c] for c in range(nT)) for g in range(len(allg))]
             for r in S]
        cc = [-sum(rows[r][c] * a0[c] for c in range(nT)) for r in S]
        t = CL.hnf_solve(B, cc)
        if t is not None:
            found = (S, t); break
    if found:
        print(f"  k = {k}: SOLVABLE rows {found[0]} => score {39033-(len(rows)-k)}")
        best = (k, found, allg, alln)
        break
    print(f"  k = {k}: no integral solution over {cnt} subsets")
pickle.dump({'gens': allg, 'gnames': alln, 'best': best, 'T': T, 'rows': rows,
             'reqs': reqs, 'a0': a0}, open(os.path.join(HERE, 'jpairs3.pkl'), 'wb'))
