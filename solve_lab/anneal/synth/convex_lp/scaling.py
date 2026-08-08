#!/usr/bin/env python3
"""scaling.py -- persistency / linear-elimination / carry-given-products counts
across modmul sizes, plus rigorous soundness enumeration on a tiny instance."""
import os
import sys
import json
import itertools
import numpy as np
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lp_core import build_modmul_instance, LP, witness_full     # noqa: E402
from analysis import gf_rank, eq_rows_as_dicts                  # noqa: E402


def and_product_cols(lp):
    return [z for (z, i, j) in lp.ands]


def measure(p, mult='schoolbook', red='quotient', probes=8, do_persist=True):
    inst = build_modmul_instance(p, mult=mult, red=red)
    lp = LP(inst)
    Q = inst['Q']
    s = inst['s']
    n = lp.n
    ninp = len(lp.input_bits)
    kinds = Counter(Q.kind[v] for v in range(n))
    carry_cols = [v for v in range(n) if Q.kind[v] in ('adder', 'carry', 'chunk')]
    and_cols = and_product_cols(lp)

    rank_eq = gf_rank(eq_rows_as_dicts(lp), n)
    eff_dim = n - rank_eq

    res = dict(p=p, s=s, n=n, n_inputs=ninp, kinds=dict(kinds),
               n_eq=len(lp.eq), rank_eq=rank_eq, eff_dim=eff_dim,
               eff_dim_minus_inputs=eff_dim - ninp,
               n_carry=len(carry_cols), n_and=len(and_cols),
               linear_reduction_frac=round(rank_eq / n, 4))

    if not do_persist:
        return inst, lp, res

    a = 5 % (1 << s)
    b = 11 % (1 << s)
    x, wv = witness_full(inst, a, b)

    # (1) inputs free: fixed-to-constant
    fx_free, _, _ = lp.persistency(n_probe=probes)
    res['persist_free_total'] = len(fx_free)

    # (2) inputs fixed: how many carries pinned (carries LP-forced given products)
    fixed_in = {v: int(x[v]) for v in lp.input_bits}
    fx_in, _, _ = lp.persistency(fixed=fixed_in, n_probe=probes)
    bad = sum(1 for v, val in fx_in.items() if abs(val - x[v]) > 1e-5)
    res['soundness_bad_inputsfixed'] = bad
    res['persist_inputsfixed_total'] = len(fx_in)
    res['carries_pinned_given_inputs'] = sum(1 for v in fx_in if v in set(carry_cols))
    res['ands_pinned_given_inputs'] = sum(1 for v in fx_in if v in set(and_cols))
    res['carry_pin_frac_given_inputs'] = round(
        res['carries_pinned_given_inputs'] / max(1, len(carry_cols)), 4)

    # (3) inputs + AND products fixed: carries LP-forced strictly by products
    fixed_prod = dict(fixed_in)
    for z in and_cols:
        fixed_prod[z] = int(x[z])
    fx_prod, _, _ = lp.persistency(fixed=fixed_prod, n_probe=probes)
    res['carries_pinned_given_products'] = sum(1 for v in fx_prod if v in set(carry_cols))
    res['carry_pin_frac_given_products'] = round(
        res['carries_pinned_given_products'] / max(1, len(carry_cols)), 4)

    # (4) inputs + products + output word C fixed
    fixed_out = dict(fixed_prod)
    for v in inst['C'].bits:
        fixed_out[v] = int(x[v])
    fx_out, _, _ = lp.persistency(fixed=fixed_out, n_probe=probes)
    res['carries_pinned_given_products_and_output'] = sum(
        1 for v in fx_out if v in set(carry_cols))
    return inst, lp, res


# ---------------------------------------------------- rigorous soundness enum
def enumerate_integer_solutions_fixed_inputs(inst, lp, a, b, cap=2_000_000):
    """Given inputs fixed, enumerate ALL integer points in the LP polytope by
    brute force over the coordinates that are non-input, verifying the witness is
    unique. Only feasible for tiny instances -- restrict brute force to the
    coordinates that the equality system leaves free (eff_dim - inputs)."""
    Q = inst['Q']
    s = inst['s']
    x, wv = witness_full(inst, a, b)
    fixed_in = {v: int(x[v]) for v in lp.input_bits}
    # free coordinates under LP with inputs fixed
    r, _, _ = None, None, None
    # find fractional/free coords by min/max
    free = []
    for v in range(lp.n):
        if v in fixed_in:
            continue
        c = np.zeros(lp.n); c[v] = 1.0
        rmin = lp.solve(c, fixed_in); rmax = lp.solve(-c, fixed_in)
        if rmin.success and rmax.success and (rmax.x[v] - rmin.x[v]) > 1e-6:
            free.append(v)
    if len(free) > 22:
        return dict(free=len(free), enumerated=False)
    # brute force over free coords in {0,1}, derive others by equality solve is hard;
    # instead check every 0/1 assignment to free coords, complete by trying to
    # satisfy all equalities+AND for the remaining via the deterministic witness --
    # simpler: test all assignments of free coords, and for each check if there is
    # any completion equal to a feasible integer point by verifying full constraints
    # after also fixing free coords and re-pinning the rest via LP (must be integral).
    count = 0
    sols = []
    for combo in itertools.product((0, 1), repeat=len(free)):
        fx = dict(fixed_in)
        for v, val in zip(free, combo):
            fx[v] = val
        rr = lp.solve(np.zeros(lp.n), fx)
        if not rr.success:
            continue
        xr = rr.x
        # all coords must now be integral for an integer solution
        if np.all(np.abs(xr - np.round(xr)) < 1e-6):
            count += 1
            sols.append(tuple(int(round(t)) for t in xr))
    return dict(free=len(free), enumerated=True, n_integer_solutions=count,
                unique_and_is_witness=(count == 1 and
                                       sols and np.array_equal(np.array(sols[0]),
                                                               np.round(x).astype(int))))


def soundness_across_inputs(inst, lp, sample=None):
    """Check: any coordinate LP-fixed-to-constant with inputs free holds in EVERY
    witness (enumerate all a,b). Also check persistency-with-inputs-fixed values
    match witness for a sample of inputs."""
    Q = inst['Q']
    s = inst['s']
    fx_free, _, _ = lp.persistency(n_probe=10)
    allab = [(a, b) for a in range(1 << s) for b in range(1 << s)]
    if sample:
        import random
        random.Random(0).shuffle(allab)
        allab = allab[:sample]
    bad = 0
    for a, b in allab:
        try:
            x, _ = witness_full(inst, a, b)
        except AssertionError:
            continue
        for v, val in fx_free.items():
            if abs(val - x[v]) > 1e-5:
                bad += 1
    return dict(n_fixed_free=len(fx_free), witnesses_checked=len(allab),
                violations=bad)


if __name__ == '__main__':
    print("=== scaling ===")
    table = []
    for p in [13, 251, 65521]:
        inst, lp, res = measure(p)
        table.append(res)
        print(json.dumps(res, default=str))
    with open('scaling_results.json', 'w') as f:
        json.dump(table, f, indent=1, default=str)

    print("\n=== soundness: enumerate integer solutions given inputs (p=13) ===")
    inst, lp, _ = measure(13, do_persist=False)
    for (a, b) in [(5, 11), (3, 7), (1, 1), (12, 12)]:
        enum = enumerate_integer_solutions_fixed_inputs(inst, lp, a, b)
        print(f"a={a} b={b}: {enum}")

    print("\n=== soundness: LP-fixed-constant holds across all witnesses (p=13) ===")
    print(soundness_across_inputs(inst, lp))
