#!/usr/bin/env python3
"""integer_gap.py -- rigorous soundness via MILP.

For a modmul instance with inputs fixed, compute the INTEGER min/max of every
variable over the feasible integer set (equalities + AND-McCormick + {0,1}).
If int-min == int-max for a variable, it is integer-determined by the inputs.
Contrast with LP min/max (convex) -> the gap is exactly what the convex
relaxation fails to pin (the carry network).
"""
import os
import sys
import json
import numpy as np
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lp_core import build_modmul_instance, LP, witness_full     # noqa: E402


def milp_minmax(lp, v, fixed, integ):
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = lp.n
    lb = np.zeros(n); ub = np.ones(n)
    for k, val in fixed.items():
        lb[k] = ub[k] = val
    cons = []
    if lp.A_eq.shape[0]:
        cons.append(LinearConstraint(lp.A_eq, lp.b_eq, lp.b_eq))
    if lp.A_ub.shape[0]:
        cons.append(LinearConstraint(lp.A_ub, -np.inf, lp.b_ub))
    c = np.zeros(n); c[v] = 1.0
    rmin = milp(c, constraints=cons, integrality=integ, bounds=Bounds(lb, ub))
    rmax = milp(-c, constraints=cons, integrality=integ, bounds=Bounds(lb, ub))
    if not (rmin.success and rmax.success):
        return None
    return rmin.x[v], -rmax.fun


def run(p, mult='schoolbook', red='quotient', a=5, b=11):
    inst = build_modmul_instance(p, mult=mult, red=red)
    lp = LP(inst)
    Q = inst['Q']
    s = inst['s']
    n = lp.n
    integ = np.ones(n)               # all binary
    a %= (1 << s); b %= (1 << s)
    x, wv = witness_full(inst, a, b)
    fixed_in = {v: int(x[v]) for v in lp.input_bits}

    int_det = 0
    int_undet = []
    lp_det = 0
    mismatch = 0
    for v in range(n):
        if v in fixed_in:
            continue
        mm = milp_minmax(lp, v, fixed_in, integ)
        if mm is None:
            continue
        imin, imax = mm
        det = abs(imax - imin) < 0.5
        if det:
            int_det += 1
            # soundness: integer-determined value equals witness
            if abs(round(imin) - x[v]) > 0.5:
                mismatch += 1
        else:
            int_undet.append(v)
        # LP min/max
        c = np.zeros(n); c[v] = 1.0
        rmn = lp.solve(c, fixed_in); rmx = lp.solve(-c, fixed_in)
        if rmn.success and rmx.success and (rmx.x[v] - rmn.x[v]) < 1e-6:
            lp_det += 1
    ninternal = n - len(fixed_in)
    kinds_undet = Counter(Q.kind[v] for v in int_undet)
    return dict(p=p, s=s, n=n, n_internal=ninternal,
                integer_determined=int_det,
                integer_undetermined=len(int_undet),
                undet_kinds=dict(kinds_undet),
                lp_determined=lp_det,
                convex_gap=int_det - lp_det,
                soundness_mismatch=mismatch,
                witness_unique=(len(int_undet) == 0))


if __name__ == '__main__':
    out = []
    for p in [13, 251]:
        r = run(p)
        out.append(r)
        print(json.dumps(r, default=str))
    with open('integer_gap.json', 'w') as f:
        json.dump(out, f, indent=1, default=str)
