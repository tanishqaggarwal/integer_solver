#!/usr/bin/env python3
"""Generic optimiser for an arbitrary target support T'.

Given T' (a set of atoms allowed to be nonzero), starting from a base state whose
nonzero atoms lie inside T':
  1. release the definers of T'  (their defined variables become free)
  2. knobs = free variables that can reach T'; rows = constraint atoms OUTSIDE T'
  3. exact integer kernel of that constraint Jacobian (forward-mode AD over Z)
  4. VERIFY each kernel vector over Z by evaluating all 39033 atoms
  5. exhaustive Hermite-normal-form search: how many of the |R(T')| rows can be
     driven to zero simultaneously by an integer point of a* + L?
  6. report the resulting score 39033 - (|R| - k)
"""
import os, pickle, sys, itertools, time
from collections import defaultdict, deque
import jengine as E, jman as J
import jcluster as CL
import jpairs as P1
import jdiag as D

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
BASE = E.load(CL.DEL)


def run(Tp, base=None, verbose=True, maxrows=26):
    Tp = sorted(set(Tp))
    Tset = set(Tp)
    val0 = list(base if base is not None else BASE)
    nz0 = [i for i in range(len(polys)) if CL.atomval(i, val0) != 0]
    assert set(nz0) <= Tset, f"base has nonzero atoms outside T': {set(nz0)-Tset}"

    definer = dict(J.definer)
    for v in [v for v, i in definer.items() if i in Tset]:
        del definer[v]
    order, cyc = E.topo(definer)
    if cyc:
        return None, "cycle after release"
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

    if forward(list(val0)) != val0:
        return None, "base not a fixed point of the released DAG"

    CONS = sorted(set(range(len(polys))) - set(definer.values()) - Tset)

    def dcolumn(z, targets):
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
                    dt = dt * val0[j] + t * der[j]
                    t *= val0[j]
                ds += dt
            der[v] = -ds // c if ds else 0
        out = {}
        for i in targets:
            s = 0
            for k, c in polys[i].items():
                t = c
                dt = 0
                for j in k:
                    dt = dt * val0[j] + t * der[j]
                    t *= val0[j]
                s += dt
            if s:
                out[i] = s
        return out

    seen = set(); q = deque()
    for i in Tp:
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

    rowsC = defaultdict(dict)
    for jx, z in enumerate(knobs):
        for i, v in dcolumn(z, CONS).items():
            rowsC[i][jx] = v
    basis = P1.int_nullspace(list(rowsC.values()), len(knobs))

    a0 = [CL.atomval(i, val0) for i in Tp]
    gens, seenv = [], set()
    for d in basis:
        v1 = list(val0)
        for k, c in enumerate(d):
            if c:
                v1[knobs[k]] += c
        forward(v1)
        nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
        if set(nz) <= Tset:
            dv = tuple(CL.atomval(i, v1) - a0[k] for k, i in enumerate(Tp))
            if any(dv) and dv not in seenv:
                seenv.add(dv); gens.append(list(dv))

    atom2eq = defaultdict(list)
    for e in CL.eqs:
        row = {}
        for c, j in e['terms']:
            row[j] = row.get(j, 0) + c
        for j in row:
            if row[j]:
                atom2eq[j].append((e['i'], row))
    R, rows = [], []
    seenR = set()
    for j in Tp:
        for (i, row) in atom2eq[j]:
            if i not in seenR:
                seenR.add(i); R.append(i)
                rows.append([row.get(x, 0) for x in Tp])
    if verbose:
        print(f"  T'={Tp}\n  knobs {len(knobs)}, kernel {len(basis)}, verified gens {len(gens)}, "
              f"|R|={len(R)}")
    if not gens or len(R) > maxrows:
        return None, f"gens={len(gens)} |R|={len(R)} (skipped)"

    nT = len(Tp)
    for k in range(len(rows), 0, -1):
        for S in itertools.combinations(range(len(rows)), k):
            B = [[sum(rows[r][c] * gens[g][c] for c in range(nT)) for g in range(len(gens))]
                 for r in S]
            cc = [-sum(rows[r][c] * a0[c] for c in range(nT)) for r in S]
            if CL.hnf_solve(B, cc) is not None:
                return (k, len(R), 39033 - (len(R) - k)), "ok"
    return (0, len(R), 39033 - len(R)), "ok"


if __name__ == '__main__':
    T = CL.all_nonzero(BASE)
    print("baseline cluster:")
    print("  ", run(T))
    for extra in ([23328], [23328, 23329], [23328, 35888], [35888], [23329]):
        Tp = sorted(set(T) | set(extra))
        print(f"\n+ atoms {extra}:")
        print("  ", run(Tp))
