#!/usr/bin/env python3
"""analysis.py -- rank / effective-dimension / TU tests on the modmul system."""
import os
import sys
import random
import numpy as np
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lp_core import build_modmul_instance, LP, witness_full   # noqa: E402


# ------------------------------------------------- exact rank over a prime field
def gf_rank(rows, ncol, mod=(1 << 61) - 1):
    """rank of an integer matrix (list of dict col->coef) over GF(mod)."""
    M = [dict((c % ncol, v % mod) for c, v in r.items()) for r in rows]
    piv_cols = {}
    rank = 0
    for r in M:
        r = dict(r)
        # reduce against existing pivots
        for c in sorted(r):
            if r.get(c, 0) == 0:
                continue
            if c in piv_cols:
                pr, pv = piv_cols[c]
                f = (r[c] * pow(pv, mod - 2, mod)) % mod
                for cc, vv in pr.items():
                    r[cc] = (r.get(cc, 0) - f * vv) % mod
                    if r[cc] == 0:
                        del r[cc]
            else:
                break
        # find leading col
        lead = None
        for c in sorted(r):
            if r.get(c, 0) != 0:
                lead = c
                break
        if lead is not None:
            piv_cols[lead] = (r, r[lead])
            rank += 1
    return rank


def eq_rows_as_dicts(lp):
    rows = []
    for lin, const in lp.eq:
        rows.append(dict(lin))
    return rows


# ------------------------------------------------- carry subsystem extraction
def carry_subsystem(lp):
    """Equality rows restricted to carry/adder columns, with product/word/input
    columns treated as fixed constants (moved to RHS). Returns (rows, carry_cols)
    where rows are dicts over carry columns only, plus the set of carry columns
    that actually appear."""
    Q = lp.Q
    carry_kinds = {'adder', 'carry', 'chunk'}
    carry_cols = set(v for v in range(lp.n) if Q.kind[v] in carry_kinds)
    rows = []
    for lin, const in lp.eq:
        row = {v: w for v, w in lin.items() if v in carry_cols}
        if row:
            rows.append(row)
    return rows, sorted(carry_cols)


def entries_hist(rows):
    c = Counter()
    for r in rows:
        for v in r.values():
            c[v] += 1
    return dict(c)


def tu_necessary(rows):
    """TU necessary condition: every entry in {-1,0,1}."""
    for r in rows:
        for v in r.values():
            if v not in (-1, 0, 1):
                return False
    return True


def tu_random_submatrix_test(rows, cols, trials=20000, kmax=6, seed=0):
    """Sample small square submatrices, compute integer determinants; TU requires
    every determinant in {-1,0,1}. Returns first violating det found (or None)."""
    rng = random.Random(seed)
    # dense-ish index
    col_index = {c: i for i, c in enumerate(cols)}
    R = [dict((col_index[c], v) for c, v in r.items()) for r in rows]
    nR, nC = len(R), len(cols)
    if nR == 0 or nC == 0:
        return ('empty', None)
    worst = 0
    for _ in range(trials):
        k = rng.randint(2, min(kmax, nR, nC))
        rs = rng.sample(range(nR), k)
        cs = rng.sample(range(nC), k)
        sub = [[R[r].get(c, 0) for c in cs] for r in rs]
        d = int(round(np.linalg.det(np.array(sub, dtype=float))))
        if abs(d) > 1:
            return ('violation', d, sub)
        worst = max(worst, abs(d))
    return ('no_violation_found', worst)


# ------------------------------------------------- fractional carry witness
def fractional_point(lp, fixed):
    """Solve LP for an arbitrary feasible point with a random objective; report
    a variable that is fractional (non-integer)."""
    rng = np.random.default_rng(1)
    r = lp.solve(rng.standard_normal(lp.n), fixed)
    frac = [(v, r.x[v]) for v in range(lp.n)
            if abs(r.x[v] - round(r.x[v])) > 1e-5]
    return r, frac


def run(p, mult='schoolbook', red='quotient'):
    inst = build_modmul_instance(p, mult=mult, red=red)
    lp = LP(inst)
    Q = inst['Q']
    s = inst['s']
    n = lp.n
    ninp = len(lp.input_bits)

    rows = eq_rows_as_dicts(lp)
    r_eq = gf_rank(rows, n)
    eff_dim = n - r_eq                      # dimension of affine hull of equalities

    csub, ccols = carry_subsystem(lp)
    hist = entries_hist(csub)
    tu_nec = tu_necessary(csub)
    tu_samp = tu_random_submatrix_test(csub, ccols)

    kinds = Counter(Q.kind[v] for v in range(n))

    out = dict(p=p, s=s, n=n, n_inputs=ninp, kinds=dict(kinds),
               n_eq=len(lp.eq), rank_eq=r_eq, eff_dim=eff_dim,
               eff_dim_minus_inputs=eff_dim - ninp,
               n_carry_cols=len(ccols), carry_entries=hist,
               tu_necessary=tu_nec, tu_sampled=tu_samp)
    return inst, lp, out


if __name__ == '__main__':
    import json
    for p in [13, 251, 65521]:
        inst, lp, out = run(p)
        print(json.dumps(out, default=str, indent=1))
        # fractional carry demonstration with inputs fixed
        s = inst['s']
        x, wv = witness_full(inst, 5 % (1 << s), 11 % (1 << s))
        fixed = {v: int(x[v]) for v in lp.input_bits}
        r, frac = fractional_point(lp, fixed)
        print(f"  inputs fixed -> LP feasible, #fractional coords = {len(frac)}; "
              f"example {frac[:3]}")
        print("---")
