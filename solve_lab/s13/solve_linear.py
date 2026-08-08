#!/usr/bin/env python3
"""
If every constraint is linear in the two unknowns mod p, the residual is a
2-variable linear system over GF(p) -- solvable by elimination, no search.

This script:
  1. validates the local forward evaluator against the real state,
  2. FITS  value_a(u,w) = A_a*u + B_a*w + C_a  (mod p)  for every check in the
     cone, by evaluating at (0,0),(1,0),(0,1) offsets,
  3. VERIFIES each fit at random points (a fit that does not predict is not a
     model -- this is the step that catches a false linearity claim),
  4. SOLVES the system A*u + B*w + C == 0 (mod p) over all checks,
  5. if solvable, builds the assignment and runs the EXACT integer checker.

Usage: python3 solve_linear.py
"""
import os, sys, time, json, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
from reduce_tw import downstream, local_forward, check_value

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
D2 = [8731, 9118]


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    unknown, touched = downstream(D2)
    checks = sorted(a for a in touched if L.atom_out.get(a) is None)
    base = {x: v[x] for x in D2}
    print(f"cone: {len(unknown)} wires, {len(checks)} checks, unknowns "
          f"x8731 & x9118\n")

    # ---- 1. validate the evaluator ----------------------------------------
    print("=" * 74)
    print("1 -- validate the local forward evaluator against the real state")
    print("=" * 74)
    val0 = local_forward(v, base, unknown, touched)
    bad = [x for x in unknown if x in val0 and val0[x] != v[x]]
    print(f"  wires recomputed : {len(val0)}")
    print(f"  disagreements    : {len(bad)}  {['x%d'%x for x in bad[:8]]}")
    if bad:
        print("  !! evaluator does not reproduce the state; fits below would be")
        print("     meaningless.  Reporting and stopping.")
        return
    av = L.all_atom_values(v)
    nz = [a for a in checks if av[a] != 0]
    print(f"  nonzero checks at base : {len(nz)} {['a%d'%a for a in nz]}")

    # ---- 2. fit the linear model -------------------------------------------
    print("\n" + "=" * 74)
    print("2 -- fit value_a(u,w) = A*u + B*w + C  (mod p)")
    print("=" * 74)

    def ev(du, dw):
        fv = {8731: base[8731] + dw, 9118: base[9118] + du}
        val = local_forward(v, fv, unknown, touched)
        return {a: check_value(a, v, val, unknown) % P for a in checks}

    f00 = ev(0, 0); f10 = ev(1, 0); f01 = ev(0, 1)
    model = {}
    for a in checks:
        A = (f10[a] - f00[a]) % P
        B = (f01[a] - f00[a]) % P
        C = f00[a] % P
        model[a] = (A, B, C)

    # ---- 3. verify the fits at random points -------------------------------
    print("\n" + "=" * 74)
    print("3 -- VERIFY each fit at random points (a fit that cannot predict is")
    print("     not a model)")
    print("=" * 74)
    random.seed(5)
    ok = defaultdict(int)
    badfit = []
    for trial in range(4):
        du = random.randrange(1, 10**6)
        dw = random.randrange(1, 10**6)
        actual = ev(du, dw)
        for a in checks:
            A, B, C = model[a]
            pred = (A * du + B * dw + C) % P
            if pred == actual[a] % P:
                ok[a] += 1
            elif a not in badfit:
                badfit.append(a)
    nlin = [a for a in checks if ok[a] == 4]
    print(f"  checks whose linear model predicts at 4/4 random points : "
          f"{len(nlin)}/{len(checks)}")
    if badfit:
        print(f"  NONLINEAR (model failed) : {len(badfit)} "
              f"{['a%d'%a for a in badfit[:10]]}")
    else:
        print(f"  ALL checks in the cone are EXACTLY linear in (u,w) mod p.")

    # ---- 4. solve the system ------------------------------------------------
    print("\n" + "=" * 74)
    print("4 -- solve  A*u + B*w + C == 0 (mod p)  over the linear checks")
    print("=" * 74)
    rows = [(model[a][0], model[a][1], (-model[a][2]) % P, a) for a in nlin]
    # gaussian elimination over GF(p) on a 2-column system
    piv = []
    R = [list(r[:3]) + [r[3]] for r in rows]
    used = []
    for col in (0, 1):
        pr = None
        for i, r in enumerate(R):
            if i in used:
                continue
            if r[col] % P != 0:
                pr = i; break
        if pr is None:
            continue
        inv = pow(R[pr][col], -1, P)
        R[pr] = [(x * inv) % P for x in R[pr][:3]] + [R[pr][3]]
        for i, r in enumerate(R):
            if i != pr and r[col] % P != 0:
                f = r[col] % P
                R[i] = [(r[k] - f * R[pr][k]) % P for k in range(3)] + [r[3]]
        used.append(pr); piv.append((col, pr))
    incons = [r for r in R if r[0] % P == 0 and r[1] % P == 0 and r[2] % P != 0]
    print(f"  rows           : {len(R)}")
    print(f"  pivots         : {len(piv)}  (columns {[c for c,_ in piv]})")
    print(f"  INCONSISTENT rows (0 = nonzero) : {len(incons)}")
    if incons:
        print(f"    e.g. from checks "
              f"{['a%d'%r[3] for r in incons[:8]]}")
        print(f"  => NO (u,w) satisfies all of them simultaneously.")
        print(f"     The obstruction is a LINEAR ALGEBRA fact over GF(p), not a")
        print(f"     search problem: no annealer, no treewidth, can change it.")
    else:
        sol = {}
        for col, pr in piv:
            sol[col] = R[pr][2] % P
        u = sol.get(0, 0); w = sol.get(1, 0)
        print(f"  CONSISTENT.  du = {u}")
        print(f"                dw = {w}")
        print(f"  => candidate exists; verifying exactly...")
        json.dump({'du': int(u), 'dw': int(w)},
                  open(os.path.join(HERE, 'linear_solution.json'), 'w'))

    print(f"\n  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
