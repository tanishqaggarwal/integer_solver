#!/usr/bin/env python3
"""Full mod-p Newton on the WHOLE constraint system, exact Jacobian by
forward-mode automatic differentiation over GF(p) (not finite differences --
three of the residual constraints are degree 2-3, where a +1 probe is wrong).

knobs = every free variable in the backward cone of the violated constraints
rows  = every constraint atom with a nonzero Jacobian entry in some knob
solve J d = -r over GF(p) by Gaussian elimination; report rank and the number of
inconsistent rows.
"""
import sys, os, time, pickle
from collections import deque
import jengine as E, jman as J, jmodp as MP, jsolve2 as S

P = MP.P
definer, order, FREE = J.definer, J.order, J.FREE
polys = E.polys
CONS = MP.CONS
EVP = MP.EVP


def cone_leaves(cs):
    seen = set(); q = deque()
    for c in cs:
        q.extend(E.varsof[c])
    lv = set()
    while q:
        x = q.popleft()
        if x in seen:
            continue
        seen.add(x)
        d = definer.get(x)
        if d is None:
            lv.add(x); continue
        for w in E.varsof[d]:
            if w != x and w not in seen:
                q.append(w)
    return lv, seen


def dcolumn(val, z):
    """Exact d/dz of every constraint atom at `val`, by forward AD mod p."""
    der = [0] * E.NV
    der[z] = 1
    for v in order:
        e = EVP.get(v)
        if e is None or v in FREE:
            continue
        ic, rest = e
        ds = 0
        for k, cc in rest:
            t = cc
            dt = 0
            for j in k:
                dt = (dt * val[j] + t * der[j]) % P
                t = t * val[j] % P
            ds += dt
        der[v] = (-ds) % P * ic % P
    out = {}
    for i in CONS:
        s = 0
        for k, c in polys[i].items():
            t = c % P
            dt = 0
            for j in k:
                dt = (dt * val[j] + t * der[j]) % P
                t = t * val[j] % P
            s += dt
        s %= P
        if s:
            out[i] = s
    return out


def gauss(rws, ncols, data, rhs):
    """Solve A d = rhs over GF(p).  data[(row,col)] -> value.
    Returns (solution dict col->val, n_inconsistent, rank)."""
    rowmap = {}
    for (i, jx), v in data.items():
        rowmap.setdefault(i, {})[jx] = v
    piv = {}
    inconsistent = 0
    for i in sorted(rws, key=lambda r: len(rowmap.get(r, {}))):
        row = dict(rowmap.get(i, {}))
        b = rhs.get(i, 0) % P
        while row:
            c = min(row)
            if c in piv:
                prow, pb = piv[c]
                f = row[c]
                for cc, vv in prow.items():
                    nv = (row.get(cc, 0) - f * vv) % P
                    if nv:
                        row[cc] = nv
                    else:
                        row.pop(cc, None)
                b = (b - f * pb) % P
            else:
                inv = pow(row[c], P - 2, P)
                piv[c] = ({cc: vv * inv % P for cc, vv in row.items()}, b * inv % P)
                row = None
                break
        if row is not None and not row and b:
            inconsistent += 1
    d = {}
    for c in sorted(piv, reverse=True):
        prow, b = piv[c]
        s = b
        for cc, vv in prow.items():
            if cc != c:
                s = (s - vv * d.get(cc, 0)) % P
        d[c] = s % P
    return d, inconsistent, len(piv)


def newton(b1, b2, iters=6):
    val, bad = S.branch(b1, b2)
    print(f"branch ({b1},{b2}) start violated: {bad}", flush=True)
    for it in range(iters):
        r0 = MP.residues(val)
        bad = sorted(i for i, x in r0.items() if x)
        print(f"\niter {it}: violated {len(bad)} {bad}", flush=True)
        if not bad:
            print("SOLVED mod p")
            return val, []
        t0 = time.time()
        lv, cone = cone_leaves(bad)
        knobs = sorted(lv)
        print(f"  cone {len(cone)} vars, knobs {len(knobs)}", flush=True)
        data = {}
        touched = set()
        for jx, z in enumerate(knobs):
            col = dcolumn(val, z)
            for i, v in col.items():
                data[(i, jx)] = v
                touched.add(i)
            if jx % 100 == 0:
                print(f"    col {jx}/{len(knobs)}  {time.time()-t0:.0f}s", flush=True)
        rws = sorted(touched)
        print(f"  Jacobian {len(rws)} x {len(knobs)}, nnz {len(data)}  ({time.time()-t0:.0f}s)", flush=True)
        missing = [i for i in bad if i not in touched]
        if missing:
            print(f"  !! violated constraints with ZERO gradient in every knob: {missing}", flush=True)
        rhs = {i: (-r0[i]) % P for i in rws}
        t1 = time.time()
        d, inc, rank = gauss(rws, len(knobs), data, rhs)
        print(f"  rank {rank} / {len(knobs)} cols, INCONSISTENT ROWS {inc}  ({time.time()-t1:.0f}s)", flush=True)
        if inc:
            print("  -> tangent system inconsistent at this point", flush=True)
            return val, bad
        for jx, z in enumerate(knobs):
            if jx in d and d[jx]:
                val[z] = (val[z] + d[jx]) % P
        MP.fwd_modp(val)
    r = MP.residues(val)
    return val, sorted(i for i, x in r.items() if x)


if __name__ == '__main__':
    b1 = int(sys.argv[1]); b2 = int(sys.argv[2])
    val, bad = newton(b1, b2)
    print("\nFINAL violated mod p:", bad)
    HERE = os.path.dirname(os.path.abspath(__file__))
    pickle.dump(val, open(os.path.join(HERE, f'jnewton2_{b1}{b2}.pkl'), 'wb'))
