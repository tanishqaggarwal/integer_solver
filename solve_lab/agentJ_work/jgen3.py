#!/usr/bin/env python3
"""Union generator search for the deliverable's 12-row / 7-atom cluster.

Two independent ways a variable can be moved without waking any atom outside T:
  (A) move it alone            -- legal iff every atom containing it lies in T
  (B) move it and let the definer DAG re-derive its forward cone
Neither dominates the other (A finds x_642/x_9413/x_17325, B finds x_8731), so take
the union, then run the exhaustive integral subset search over the resulting lattice.

score = 39033 - (12 - k), k = number of the 12 rows driven to zero.  k=5 -> 39026.
"""
import os, pickle, sys, itertools
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL
import jgen2 as G2

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys


def collect(val0, T, cands):
    Tset = set(T)
    a0 = [CL.atomval(i, val0) for i in T]
    gens, gnames = [], []
    seen = set()
    for v in sorted(cands):
        for mode in ('alone', 'dag'):
            v1 = list(val0); v1[v] += 1
            if mode == 'dag':
                G2.fwd_from(v1, v)
            nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
            if set(nz) <= Tset:
                d1 = tuple(CL.atomval(i, v1) - a0[k] for k, i in enumerate(T))
                if any(d1) and d1 not in seen:
                    v2 = list(val0); v2[v] += 2
                    if mode == 'dag':
                        G2.fwd_from(v2, v)
                    d2 = [CL.atomval(i, v2) - a0[k] for k, i in enumerate(T)]
                    lin = all(d2[k] == 2 * d1[k] for k in range(len(T)))
                    seen.add(d1)
                    gens.append(list(d1))
                    gnames.append(f"x_{v}[{mode}]{'' if lin else '(NONLIN)'}")
    return gens, gnames, a0


if __name__ == '__main__':
    val0 = E.load(CL.DEL)
    T = CL.all_nonzero(val0)
    Tset = set(T)
    print("T =", T)

    # broad candidate set: variables of T's atoms, their backward cone, and any
    # variable sharing an atom with those (2-hop in the variable/atom graph)
    base = set()
    for i in T:
        base |= E.varsof[i]
    seen = set(); q = deque(base)
    while q:
        x = q.popleft()
        if x in seen: continue
        seen.add(x)
        d = J.definer.get(x)
        if d is None: continue
        for w in E.varsof[d]:
            if w != x and w not in seen: q.append(w)
    cands = set(seen)
    for v in list(seen):
        for i in CL.occ[v]:
            cands |= E.varsof[i]
    cands |= {v for v in range(E.NV) if set(CL.occ[v]) <= Tset}
    print("candidates:", len(cands))

    gens, gnames, a0 = collect(val0, T, cands)
    print(f"\nUNION generators ({len(gens)}):")
    for n_, g in zip(gnames, gens):
        print("   ", n_, [x if abs(x) < 10**10 else f"~{len(str(abs(x)))}d" for x in g])

    rows, reqs = [], []
    for e in CL.eqs:
        row = {}
        for c, j in e['terms']:
            row[j] = row.get(j, 0) + c
        if any(j in Tset and row[j] for j in row):
            reqs.append(e['i']); rows.append([row.get(j, 0) for j in T])
    print("rows:", len(rows), reqs)
    cur = [sum(rows[r][k] * a0[k] for k in range(len(T))) for r in range(len(rows))]
    print("currently zero rows:", sum(1 for x in cur if x == 0))

    nT = len(T)
    best = None
    for k in range(len(rows), 4, -1):
        found = None; cnt = 0
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
            best = (k, found); break
        print(f"  k = {k}: no integral solution over {cnt} subsets")
    pickle.dump({'T': T, 'rows': rows, 'reqs': reqs, 'a0': a0, 'gens': gens,
                 'gnames': gnames, 'best': best},
                open(os.path.join(HERE, 'jgen3.pkl'), 'wb'))
