#!/usr/bin/env python3
"""Compensated moves: enlarge the reachable lattice L on the deliverable's cluster.

Single-variable moves (alone / with forward re-derivation) gave 10 generators and the
exhaustive HNF search capped at k=5.  The untested class is JOINT moves: change u and
v (and more) together so their effects on every atom OUTSIDE the support T cancel,
while their effect on the 7 cluster values does not.

Exactness discipline: atoms have degree <= 2, so an atom is EXACTLY affine along a
move supported on a set S provided no two members of S co-occur in a quadratic
monomial of that atom.  We enforce that, so:
  * gradients are exact integers,
  * an integer kernel vector is valid at ANY integer multiple (no second-order leak),
  * and we still verify every candidate by evaluating all 39033 atoms over Z.
"""
import os, pickle, sys, itertools
from collections import defaultdict, deque
from fractions import Fraction
import jengine as E, jman as J
import jcluster as CL

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
occ = CL.occ


def quad_pairs(i):
    """unordered variable pairs occurring in a quadratic monomial of atom i"""
    out = set()
    for k in polys[i]:
        if len(k) == 2:
            out.add((min(k), max(k)))
    return out


def grad(i, v, val):
    """exact d(atom i)/d(x_v) at val"""
    g = 0
    for k, c in polys[i].items():
        if v not in k:
            continue
        if len(k) == 1:
            g += c
        else:
            a, b = k
            if a == b:
                g += 2 * c * val[a]
            else:
                other = b if a == v else a
                g += c * val[other]
    return g


def int_nullspace(rows, n):
    """Integer nullspace basis of the matrix given by `rows` (list of dict col->val)."""
    # fraction-free row reduction, tracking column ops implicitly via rational solve
    A = [[Fraction(r.get(c, 0)) for c in range(n)] for r in rows]
    m = len(A)
    piv = []
    r = 0
    for c in range(n):
        p = None
        for i in range(r, m):
            if A[i][c]:
                p = i; break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [A[i][k] - f * A[r][k] for k in range(n)]
        piv.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        vec = [Fraction(0)] * n
        vec[fc] = Fraction(1)
        for i, c in enumerate(piv):
            vec[c] = -A[i][fc]
        den = 1
        for x in vec:
            den = den * x.denominator // __import__('math').gcd(den, x.denominator)
        basis.append([int(x * den) for x in vec])
    return basis


if __name__ == '__main__':
    val = E.load(CL.DEL)
    T = CL.all_nonzero(val)
    Tset = set(T)
    a0 = [CL.atomval(i, val) for i in T]
    print("T =", T)

    # variable cone: variables of T's atoms, plus one hop through shared atoms
    S = set()
    for i in T:
        S |= E.varsof[i]
    hop = set(S)
    for v in list(S):
        for i in occ[v]:
            hop |= E.varsof[i]
    S = hop
    print("candidate move variables (1 hop):", len(S))

    # atoms outside T touched by S
    OUT = sorted({i for v in S for i in occ[v]} - Tset)
    print("outside-T atoms touched:", len(OUT))

    # enforce exact affineness: drop variables that co-occur quadratically with
    # another candidate inside any touched atom (greedy independent set)
    conflict = defaultdict(set)
    for i in OUT + T:
        for (a, b) in quad_pairs(i):
            if a in S and b in S:
                conflict[a].add(b); conflict[b].add(a)
    order = sorted(S, key=lambda v: len(conflict[v]))
    keep = []
    used = set()
    for v in order:
        if v in used:
            continue
        keep.append(v); used.add(v); used |= conflict[v]
    S = sorted(keep)
    print("after quadratic-independence filter:", len(S))

    OUT = sorted({i for v in S for i in occ[v]} - Tset)
    print("outside-T atoms in play:", len(OUT))

    idx = {v: k for k, v in enumerate(S)}
    rows = []
    for i in OUT:
        row = {}
        for v in E.varsof[i]:
            if v in idx:
                g = grad(i, v, val)
                if g:
                    row[idx[v]] = g
        if row:
            rows.append(row)
    print("gradient matrix:", len(rows), "x", len(S))

    basis = int_nullspace(rows, len(S))
    print("integer kernel dimension:", len(basis))

    gens, gnames = [], []
    seen = set()
    for bi, d in enumerate(basis):
        v1 = list(val)
        for k, c in enumerate(d):
            if c:
                v1[S[k]] += c
        nz = [i for i in range(len(polys)) if CL.atomval(i, v1) != 0]
        if set(nz) <= Tset:
            dv = tuple(CL.atomval(i, v1) - a0[k] for k, i in enumerate(T))
            if any(dv) and dv not in seen:
                seen.add(dv)
                gens.append(list(dv))
                sup = [S[k] for k, c in enumerate(d) if c]
                gnames.append(f"ker{bi}[{len(sup)} vars]")
        else:
            pass
    print(f"\nVERIFIED compensated generators: {len(gens)}")
    for n_, g in zip(gnames, gens):
        print("   ", n_, [x if abs(x) < 10**10 else f"~{len(str(abs(x)))}d" for x in g])
    pickle.dump({'gens': gens, 'gnames': gnames, 'S': S},
                open(os.path.join(HERE, 'jpairs.pkl'), 'wb'))
