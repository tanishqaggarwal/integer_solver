#!/usr/bin/env python3
"""Widen the generator set for the deliverable's cluster.

A variable is usable if, after perturbing it and letting the definer DAG
re-derive its forward cone, the set of nonzero atoms is still inside T.
That admits compensating chains (x_9118 / x_8731 were claimed usable this way),
which the strictly-safe test in jcluster.py rejects.
Everything is checked empirically -- the test "nonzero atoms subset of T" is exact.
"""
import os, pickle, sys, itertools, random
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
definer, order, FREE = J.definer, J.order, J.FREE
pos = {v: k for k, v in enumerate(order)}
uses = defaultdict(list)
for v, i in definer.items():
    for w in E.varsof[i]:
        if w != v:
            uses[w].append(v)
EV = J.ev


def fwd_from(val, seed):
    dirty = set(); q = deque([seed])
    while q:
        x = q.popleft()
        for w in uses[x]:
            if w not in dirty and w not in FREE:
                dirty.add(w); q.append(w)
    for v in sorted(dirty, key=lambda z: pos[z]):
        c, rest = EV[v]
        s = 0
        for k, cc in rest:
            t = cc
            for j in k:
                t *= val[j]
            s += t
        val[v] = (-s) // c
    return val


if __name__ == '__main__':
    val0 = E.load(CL.DEL)
    T = CL.all_nonzero(val0)
    Tset = set(T)
    a0 = [CL.atomval(i, val0) for i in T]
    print("T =", T)

    # candidate variables: backward cone of T through the definer DAG
    seen = set(); q = deque()
    for i in T:
        q.extend(E.varsof[i])
    cand = set()
    while q:
        x = q.popleft()
        if x in seen: continue
        seen.add(x); cand.add(x)
        d = definer.get(x)
        if d is None: continue
        for w in E.varsof[d]:
            if w != x and w not in seen: q.append(w)
    cand |= {v for v in range(E.NV) if set(CL.occ[v]) <= Tset}
    print("candidate variables:", len(cand))

    gens, gnames = [], []
    for v in sorted(cand):
        v1 = list(val0); v1[v] += 1
        fwd_from(v1, v)
        nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
        if set(nz) <= Tset:
            d1 = [CL.atomval(i, v1) - a0[k] for k, i in enumerate(T)]
            # linearity check at +2
            v2 = list(val0); v2[v] += 2
            fwd_from(v2, v)
            d2 = [CL.atomval(i, v2) - a0[k] for k, i in enumerate(T)]
            lin = all(d2[k] == 2 * d1[k] for k in range(len(T)))
            if any(d1):
                gens.append(d1); gnames.append(f"x_{v}{'' if lin else ' (NONLINEAR)'}")
    print(f"\nusable generators ({len(gens)}):")
    for n_, g in zip(gnames, gens):
        print("   ", n_, [x if abs(x) < 10**12 else f"~{len(str(abs(x)))}d" for x in g])

    rows, reqs = CL.__dict__.get('_rows'), None
    # rebuild rows
    rows = []; reqs = []
    for e in CL.eqs:
        row = {}
        for c, j in e['terms']:
            row[j] = row.get(j, 0) + c
        if any(j in Tset and row[j] for j in row):
            reqs.append(e['i']); rows.append([row.get(j, 0) for j in T])
    print(f"rows {len(rows)}")

    nT = len(T)
    best = None
    for k in range(len(rows), 4, -1):
        found = None
        cnt = 0
        for S in itertools.combinations(range(len(rows)), k):
            cnt += 1
            B = [[sum(rows[r][c] * gens[g][c] for c in range(nT)) for g in range(len(gens))]
                 for r in S]
            cc = [-sum(rows[r][c] * a0[c] for c in range(nT)) for r in S]
            t = CL.hnf_solve(B, cc)
            if t is not None:
                found = (S, t); break
        if found:
            print(f"  k = {k}: SOLVABLE rows {found[0]} => score {39033-(len(rows)-k)}")
            best = (k, found, gens, gnames)
            break
        print(f"  k = {k}: no integral solution over {cnt} subsets")
    pickle.dump({'T': T, 'rows': rows, 'reqs': reqs, 'a0': a0, 'gens': gens,
                 'gnames': gnames, 'best': best},
                open(os.path.join(HERE, 'jgen2.pkl'), 'wb'))
