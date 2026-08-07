#!/usr/bin/env python3
"""Mod-p Newton on the full constraint system.

Knobs: every free variable that moves any currently-violated constraint.
Rows:  every constraint that any knob moves.
Build the exact Jacobian by probing (the map is polynomial; the probe gives the
finite difference, which for the linear part is exact), solve J d = -r mod p by
Gaussian elimination, apply, iterate.
"""
import sys, random, time
from collections import deque
import jengine as E, jman as J, jmodp as MP, jsolve2 as S

P = MP.P
definer = J.definer


def cone_leaves(cs):
    seen = set(); q = deque()
    for c in cs:
        q.extend(E.varsof[c])
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
    return lv


def jac(val, knobs, rows=None):
    """returns (rows, cols, matrix as dict[(r,c)] ) using +1 finite differences."""
    r0 = MP.residues(val)
    cols = list(knobs)
    data = {}
    touched = set()
    for jx, z in enumerate(cols):
        v2 = list(val); v2[z] = (v2[z] + 1) % P
        MP.fwd_modp(v2)
        r2 = MP.residues(v2)
        for i in MP.CONS:
            d = (r2[i] - r0[i]) % P
            if d:
                data[(i, jx)] = d
                touched.add(i)
    rws = sorted(touched if rows is None else set(rows) | touched)
    return rws, cols, data, r0


def gauss_solve(rws, cols, data, rhs):
    """Solve A d = rhs mod p.  Returns (d, n_inconsistent)."""
    n = len(cols)
    rowsd = []
    for ri, i in enumerate(rws):
        row = {}
        for jx in range(n):
            v = data.get((i, jx))
            if v: row[jx] = v
        rowsd.append((row, rhs.get(i, 0) % P))
    piv = {}
    order_rows = sorted(range(len(rowsd)), key=lambda k: len(rowsd[k][0]))
    reduced = []
    for k in order_rows:
        row, b = rowsd[k]
        row = dict(row)
        for c, (prow, pb) in piv.items():
            if c in row:
                f = row[c]
                for cc, vv in prow.items():
                    row[cc] = (row.get(cc, 0) - f * vv) % P
                    if row[cc] == 0: del row[cc]
                b = (b - f * pb) % P
        if not row:
            reduced.append(('inconsistent' if b else 'ok', k))
            continue
        c = min(row)
        inv = pow(row[c], P - 2, P)
        prow = {cc: vv * inv % P for cc, vv in row.items()}
        piv[c] = (prow, b * inv % P)
        reduced.append(('pivot', k))
    inc = sum(1 for t, _ in reduced if t == 'inconsistent')
    # back substitution: free columns = 0
    d = {}
    for c in sorted(piv, reverse=True):
        prow, b = piv[c]
        s = b
        for cc, vv in prow.items():
            if cc != c:
                s = (s - vv * d.get(cc, 0)) % P
        d[c] = s % P
    return d, inc, len(piv)


if __name__ == '__main__':
    b1 = int(sys.argv[1]); b2 = int(sys.argv[2])
    val, bad = S.branch(b1, b2)
    print("start violated:", bad)
    for it in range(6):
        r0 = MP.residues(val)
        bad = [i for i, x in r0.items() if x]
        print(f"iter {it}: violated {len(bad)} {bad}")
        if not bad:
            break
        t0 = time.time()
        knobs = sorted(cone_leaves(bad))
        print(f"  knobs {len(knobs)}")
        rws, cols, data, r0 = jac(val, knobs)
        print(f"  jacobian {len(rws)} x {len(cols)} nnz {len(data)}  ({time.time()-t0:.0f}s)")
        rhs = {i: (-r0[i]) % P for i in rws}
        d, inc, rank = gauss_solve(rws, cols, data, rhs)
        print(f"  rank {rank}, inconsistent rows {inc}")
        if inc:
            print("  -> tangent system inconsistent; stopping")
            break
        for jx, z in enumerate(cols):
            if jx in d:
                val[z] = (val[z] + d[jx]) % P
        MP.fwd_modp(val)
    r = MP.residues(val)
    print("final violated:", sorted(i for i, x in r.items() if x))
    import pickle, os
    pickle.dump(val, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       f'jnewton_{b1}{b2}.pkl'), 'wb'))
