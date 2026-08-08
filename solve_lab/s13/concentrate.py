#!/usr/bin/env python3
"""
CONCENTRATE THE HARDNESS.

Established this session:
  * given a frame, every residual check is EXACTLY LINEAR mod p in the
    continuous knobs (solve_linear3.py, verified by prediction);
  * so the continuous part of the problem is linear algebra -- easy, and NOT
    where the difficulty lives.

Therefore the difficulty must live in the DISCRETE part: the message bits, which
enter through load pins of the form   bit*(x_B - HUGE) - s*x_C.  If each set bit
adds a fixed residue, then every downstream quantity -- including the obstruction
-- is AFFINE in the bit vector over GF(p):

        INV(b)  ==  INV(0) + sum_i b_i * v_i   (mod p)

If that holds, the whole instance concentrates into

        find b in {0,1}^n  with  sum_i b_i * v_i == t  (mod p)

a MODULAR SUBSET-SUM: n binary unknowns, a handful of 256-bit congruences.  That
is small, dense, genuinely hard, and exactly the shape a QUBO/annealer wants --
the hardness moved INTO the block instead of into the coupling.

This script MEASURES whether that affinity holds, rather than assuming it:
  1. find the boolean free inputs that actually move the residual,
  2. measure each one's individual effect v_i,
  3. TEST additivity on random subsets (the step that decides everything),
  4. if affine, emit the subset-sum instance and size its QUBO.

Usage: python3 concentrate.py
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
NA = len(L.avars)


def residual_signature(v):
    """The obstruction: values mod p of every CHECK atom that is nonzero."""
    av = L.all_atom_values(v)
    sig = {}
    for a in range(NA):
        if L.atom_out.get(a) is None:
            r = av[a] % P
            if r:
                sig[a] = r
    return sig


def forward_from_bits(base, bits_on, boolvars):
    """Set the boolean inputs, clean-forward the whole circuit, return the state."""
    v = list(base)
    for b in boolvars:
        v[b] = 1 if b in bits_on else 0
    ad.fwd(v, rounds=6)
    return v


def main():
    t0 = time.time()
    v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    free = [t for t in range(L.NVARS) if t not in L.definer]

    # ---- 1. the boolean free inputs ---------------------------------------
    boolvars = [x for x in free if v0[x] in (0, 1)]
    on_now = [x for x in boolvars if v0[x] == 1]
    print(f"boolean free inputs : {len(boolvars):,}   currently set: {len(on_now)} "
          f"{['x%d'%x for x in on_now[:8]]}")

    # baseline: clean forward from the current bits
    t1 = time.time()
    vb = forward_from_bits(v0, set(on_now), boolvars)
    sig0 = residual_signature(vb)
    print(f"clean forward: {time.time()-t1:.2f}s; "
          f"nonzero checks at baseline: {len(sig0)}")
    if not sig0:
        print("  baseline already has an empty residual -- nothing to concentrate")
        return
    targets = sorted(sig0)[:6]
    print(f"tracking obstruction on checks {['a%d'%a for a in targets]}")

    # ---- 2. per-bit effect -------------------------------------------------
    print("\n" + "=" * 74)
    print("MEASURE the individual effect of each bit (sampled)")
    print("=" * 74)
    random.seed(2)
    cand = [b for b in boolvars if b not in set(on_now)]
    sample = random.sample(cand, min(40, len(cand)))
    eff = {}
    live = []
    for i, b in enumerate(sample):
        vv = forward_from_bits(v0, set(on_now) | {b}, boolvars)
        sg = residual_signature(vv)
        d = {a: (sg.get(a, 0) - sig0.get(a, 0)) % P for a in targets}
        eff[b] = d
        if any(d.values()):
            live.append(b)
    print(f"  bits sampled            : {len(sample)}")
    print(f"  bits that MOVE the obstruction : {len(live)} "
          f"{['x%d'%b for b in live[:10]]}")
    if not live:
        print("  none of the sampled bits moves it -- affinity test not applicable")
        print("  (this matches the lab's finding that most booleans are inert)")
        return

    # ---- 3. additivity test ------------------------------------------------
    print("\n" + "=" * 74)
    print("TEST ADDITIVITY:  INV(b_i + b_j) - INV(0) == v_i + v_j  (mod p) ?")
    print("=" * 74)
    ok = bad = 0
    for _ in range(min(8, len(live) * (len(live) - 1) // 2 or 1)):
        if len(live) < 2:
            break
        bi, bj = random.sample(live, 2)
        vv = forward_from_bits(v0, set(on_now) | {bi, bj}, boolvars)
        sg = residual_signature(vv)
        good = True
        for a in targets:
            actual = (sg.get(a, 0) - sig0.get(a, 0)) % P
            pred = (eff[bi][a] + eff[bj][a]) % P
            if actual != pred:
                good = False
        (ok, bad) = (ok + 1, bad) if good else (ok, bad + 1)
        print(f"  pair (x{bi}, x{bj}): {'ADDITIVE' if good else 'NOT additive'}")
    print(f"\n  additive pairs: {ok}   non-additive: {bad}")

    if bad == 0 and ok > 0:
        print("""
  => the obstruction IS affine in the bits over GF(p).
     The instance concentrates into a MODULAR SUBSET-SUM:
         find b in {0,1}^n with sum_i b_i * v_i == t (mod p)
     n binary unknowns, one 256-bit congruence per target.""")
        inst = {'p': str(P),
                'targets': {f'a{a}': str((-sig0[a]) % P) for a in targets},
                'coeffs': {f'x{b}': {f'a{a}': str(eff[b][a]) for a in targets}
                           for b in live}}
        out = os.path.join(HERE, 'subsetsum.json')
        json.dump(inst, open(out, 'w'), indent=1)
        print(f"     instance written -> {out}")
        n = len(live)
        print(f"\n  QUBO SIZE of the concentrated core:")
        print(f"     binary unknowns (bits)        : {n}")
        print(f"     congruences                   : {len(targets)} x 256 bits")
        print(f"     carry/limb auxiliaries        : "
              f"~{len(targets)} x 16 limbs x ~50 = {len(targets)*16*50:,}")
        print(f"     TOTAL                         : ~{n + len(targets)*16*50:,} "
              f"binary  -- one block, no coupling")
    else:
        print("""
  => NOT affine.  The bits interact, so the concentrated object is a
     higher-degree boolean system, not a plain subset-sum.  Report only.""")

    print(f"\n  {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
