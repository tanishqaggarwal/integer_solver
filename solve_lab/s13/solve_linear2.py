#!/usr/bin/env python3
"""
Redo of solve_linear.py with the BLOCK-PRESERVING evaluator (fwd_frame.py).

The first attempt used a naive forward pass that silently "repaired" the two
deliberately-broken gates of the 39,026 witness, changing 109 downstream wires.
Its linearity verdict (60/60 linear) was computed on wrong values and is
WITHDRAWN.  Everything below uses an evaluator whose identity test is exact.

Steps: fit -> VERIFY the fit at random points -> solve -> exact integer check.

Usage: python3 solve_linear2.py
"""
import os, sys, time, json, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
from reduce_tw import downstream
import fwd_frame as F

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
D2 = [8731, 9118]


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    unknown, touched = downstream(D2)
    checks = sorted(a for a in touched if L.atom_out.get(a) is None)
    frozen, bad = F.self_test(v, unknown, touched)
    print(f"cone {len(unknown)} wires / {len(checks)} checks; "
          f"frozen {len(frozen)}; identity test disagreements {len(bad)}")
    assert not bad, "evaluator must reproduce the state"
    base = {8731: v[8731], 9118: v[9118]}

    def ev(du, dw):
        ov = {9118: base[9118] + du, 8731: base[8731] + dw}
        val = F.evaluate(v, ov, unknown, touched, frozen)
        return {a: F.check_value(a, v, val) for a in checks}

    # ---------------- fit ----------------------------------------------------
    print("\n" + "=" * 74)
    print("FIT  value_a(u,w) = A*u + B*w + C   (mod p)")
    print("=" * 74)
    f00, f10, f01 = ev(0, 0), ev(1, 0), ev(0, 1)
    model = {a: ((f10[a] - f00[a]) % P, (f01[a] - f00[a]) % P, f00[a] % P)
             for a in checks}
    nz0 = [a for a in checks if f00[a] != 0]
    print(f"  checks nonzero at base : {len(nz0)} {['a%d'%a for a in nz0]}")

    # ---------------- verify the fit ----------------------------------------
    print("\n" + "=" * 74)
    print("VERIFY the fit at random points (this is what caught the last bug)")
    print("=" * 74)
    random.seed(9)
    good = defaultdict(int)
    trials = 4
    for _ in range(trials):
        du, dw = random.randrange(1, 10**9), random.randrange(1, 10**9)
        act = ev(du, dw)
        for a in checks:
            A, B, C = model[a]
            if (A * du + B * dw + C) % P == act[a] % P:
                good[a] += 1
    lin = [a for a in checks if good[a] == trials]
    non = [a for a in checks if good[a] != trials]
    print(f"  linear mod p (4/4 predictions correct) : {len(lin)}/{len(checks)}")
    print(f"  NOT linear                             : {len(non)} "
          f"{['a%d'%a for a in non[:12]]}")

    # ---------------- solve over the linear rows -----------------------------
    print("\n" + "=" * 74)
    print("SOLVE  A*u + B*w + C == 0 (mod p) over the linear checks")
    print("=" * 74)
    rows = [[model[a][0] % P, model[a][1] % P, (-model[a][2]) % P, a]
            for a in lin]
    R = [r[:] for r in rows]
    piv, used = [], set()
    for col in (0, 1):
        pr = None
        for i, r in enumerate(R):
            if i in used or r[col] % P == 0:
                continue
            pr = i; break
        if pr is None:
            continue
        inv = pow(R[pr][col], -1, P)
        R[pr] = [(x * inv) % P for x in R[pr][:3]] + [R[pr][3]]
        for i, r in enumerate(R):
            if i != pr and r[col] % P != 0:
                f = r[col] % P
                R[i] = [(r[k] - f * R[pr][k]) % P for k in range(3)] + [r[3]]
        used.add(pr); piv.append((col, pr))
    incons = [r for r in R if r[0] % P == 0 and r[1] % P == 0 and r[2] % P != 0]
    print(f"  linear rows : {len(R)}   pivots : {len(piv)}   "
          f"inconsistent rows : {len(incons)}")

    if incons:
        print(f"    inconsistent from checks "
              f"{['a%d'%r[3] for r in incons[:10]]}")
        print("  => within this frame no (u,w) satisfies all linear checks.")
    else:
        sol = {c: R[pr][2] % P for c, pr in piv}
        du, dw = sol.get(0, 0), sol.get(1, 0)
        print(f"  CONSISTENT: du = {du}\n              dw = {dw}")
        # exact verification of the resulting full assignment
        ov = {9118: base[9118] + du, 8731: base[8731] + dw}
        val = F.evaluate(v, ov, unknown, touched, frozen)
        newv = list(v)
        for x, t in val.items():
            newv[x] = t
        out = os.path.join(HERE, 'cand_linear.json')
        json.dump({f'x_{i}': int(newv[i]) for i in range(len(newv)) if newv[i]},
                  open(out, 'w'))
        print(f"  candidate written -> {out}")
        print(f"  verify with: python3 {LAB}/checker.py {out}")

    print(f"\n  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
