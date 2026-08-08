#!/usr/bin/env python3
"""
CONCENTRATION, second axis: push the hardness into the QUADRATIC CORE.

What this session established:
  * per frame, every residual check is LINEAR mod p in the continuous knobs
    -> the continuous part is elimination, not search;
  * the message bits are NOT additive (concentrate.py) -> they are not a
    subset-sum axis either.

So the difficulty is neither the linear part nor a plain boolean sum.  What is
left is the ALGEBRAIC CORE that REDUCED_PROBLEM.md section 3 isolated:

    S = A*u^2 - w^2        T = B*u - w*c          require S == T == 0 (mod p)
    u=x29322  w=x3558  A=x33469  B=x27713  c=x1326

Eliminating w gives  u^2 * (A*c^2 - B^2) == 0 (mod p), so either u == w == 0
(the degenerate branch the witness sits in) or

        A * c^2  ==  B^2   (mod p)          <-- ONE quadratic condition

That single congruence is the natural place to concentrate: it is small (five
256-bit unknowns), it is genuinely hard (a quadratic-residue / square-root
condition in GF(p)), and everything else around it is linear.

This script MEASURES, at the current state:
  1. the core quantities and whether the core identity really controls the checks,
  2. how A, B, c depend on the knobs (affine mod p?), since that decides whether
     the concentrated QUBO is over a few residues or over the whole circuit,
  3. the QUBO size of the concentrated core.

Usage: python3 concentrate2.py
"""
import os, sys, time, json, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, os.path.join(HERE, '..', 's10'))
sys.path.insert(0, HERE)
import lib as L
import ad

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
CORE = {'u': 29322, 'w': 3558, 'A': 33469, 'B': 27713, 'c': 1326}


def leg(a):
    a %= P
    return 0 if a == 0 else (1 if pow(a, (P - 1) // 2, P) == 1 else -1)


def main():
    t0 = time.time()
    v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    v = list(v0); ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    nz = [a for a in range(len(L.avars))
          if L.atom_out.get(a) is None and av[a] != 0]
    print(f"clean frame: score {L.NEQ - len(L.failing_eqs(av))}, "
          f"{len(nz)} obstruction checks {['a%d'%a for a in nz]}")

    # ---- 1. the core quantities -------------------------------------------
    print("\n" + "=" * 74)
    print("1 -- the algebraic core at this state")
    print("=" * 74)
    val = {k: v[i] % P for k, i in CORE.items()}
    for k, i in CORE.items():
        s = str(val[k])
        print(f"  {k} = x{i:<6} = {s[:40]}{'...' if len(s) > 40 else ''}")
    A, B, c, u, w = (val[k] for k in ('A', 'B', 'c', 'u', 'w'))
    S = (A * u * u - w * w) % P
    T = (B * u - w * c) % P
    disc = (A * c * c - B * B) % P
    print(f"\n  S = A u^2 - w^2 = {'0' if S == 0 else 'NONZERO'}")
    print(f"  T = B u - w c   = {'0' if T == 0 else 'NONZERO'}")
    print(f"  u = {'0' if u == 0 else 'nonzero'},  w = {'0' if w == 0 else 'nonzero'}"
          f"   -> the witness sits in the DEGENERATE branch u=w=0"
          if u == 0 and w == 0 else "")
    print(f"\n  the non-degenerate branch needs  A*c^2 == B^2 (mod p):")
    print(f"    A*c^2 - B^2 = {'0  (SATISFIED)' if disc == 0 else 'NONZERO'}")
    print(f"    Legendre (A|p) = {leg(A):+d}   "
          f"({'A is a QR: a square root exists' if leg(A) == 1 else 'A is a NON-residue: this branch is closed'})")
    if c % P:
        rhs = (B * pow(c, -1, P)) % P
        print(f"    A == (B/c)^2 ? {A == rhs * rhs % P}")

    # ---- 2. how do A, B, c depend on the knobs? ----------------------------
    print("\n" + "=" * 74)
    print("2 -- are the core quantities AFFINE in the knobs (mod p)?")
    print("=" * 74)
    free = [t for t in range(L.NVARS) if t not in L.definer]
    bigfree = [x for x in free if 200 <= abs(v0[x]).bit_length() <= 400]
    print(f"  large free inputs (advice-scale): {len(bigfree)} "
          f"{['x%d'%x for x in bigfree[:14]]}")

    def core_at(delta):
        vv = list(v0)
        for x, d in delta.items():
            vv[x] += d
        ad.fwd(vv, rounds=6)
        return {k: vv[i] % P for k, i in CORE.items()}

    base = core_at({})
    eff = {}
    for x in bigfree[:10]:
        e1 = core_at({x: 1})
        eff[x] = {k: (e1[k] - base[k]) % P for k in CORE}
    movers = [x for x in eff if any(eff[x].values())]
    print(f"  knobs that move the core: {len(movers)} {['x%d'%x for x in movers]}")

    add_ok = add_bad = 0
    if len(movers) >= 2:
        random.seed(1)
        for _ in range(4):
            bi, bj = random.sample(movers, 2)
            a2 = core_at({bi: 1, bj: 1})
            good = all((a2[k] - base[k]) % P ==
                       (eff[bi][k] + eff[bj][k]) % P for k in CORE)
            (add_ok, add_bad) = (add_ok + 1, add_bad) if good else (add_ok, add_bad + 1)
        print(f"  additivity over knob pairs: {add_ok} additive, {add_bad} not")
    affine = add_bad == 0 and add_ok > 0

    # ---- 3. QUBO size of the concentrated core ----------------------------
    print("\n" + "=" * 74)
    print("3 -- QUBO SIZE of the concentrated core")
    print("=" * 74)
    MUL = 39750          # verified 256x256 modular multiply, 16 blocks
    print(f"  condition            : A*c^2 == B^2 (mod p)")
    print(f"  modular multiplies   : c^2, A*c^2, B^2  = 3")
    print(f"  per multiply         : {MUL:,} binary (16 blocks of ~2,484)")
    print(f"  core QUBO            : ~{3*MUL:,} binary in ~48 blocks")
    if affine:
        print(f"  unknowns             : the {len(movers)} knobs are affine in the")
        print(f"                         core, so A,B,c are affine forms -> the")
        print(f"                         QUBO is over {len(movers)} x 256 = "
              f"{len(movers)*256:,} binary of true freedom")
    else:
        print(f"  unknowns             : core NOT affine in the knobs, so A,B,c")
        print(f"                         must be computed through the circuit;")
        print(f"                         add the cone's arithmetic to the block")
    print(f"\n  => this is the concentrated object: ONE dense hard block of")
    print(f"     ~{3*MUL:,} binary variables carrying a genuine GF(p) quadratic")
    print(f"     condition, with everything around it linear (eliminable).")
    print(f"     Contrast the earlier decomposition, whose blocks were all")
    print(f"     deterministic circuit evaluations with the hardness in the")
    print(f"     coupling instead.")

    print(f"\n  {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
