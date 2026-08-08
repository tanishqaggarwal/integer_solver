#!/usr/bin/env python3
"""
The EXTENDED core: the 2 residual conditions plus the collateral that makes them
hard, sized for QUBO decomposition.

Decoded core (core_print.py), with x26064 = p and x1329, x10903 free handles:

    a35759 = 0   <=>  5113045 * x7075 * x9118 = x1329  * p   <=>  p | x7075*x9118
    a35760 = 0   <=>            x7075 * x8731 = -x10903 * p  <=>  p | x7075*x8731
    x7075 = 1 - x2081*x4287        (x2081, x4287 boolean)

So the residual is TWO divisibilities.  Either kill the selector (x2081=x4287=1
=> x7075=0) or make x9118 and x8731 multiples of p.  Both are free inputs, so the
conditions are trivially reachable IN ISOLATION -- the difficulty is COLLATERAL:
x9118, x8731, x2081, x4287 occur in other checks.

This script measures that collateral exactly:
  * every atom containing one of the 4 driver free inputs,
  * which are checks (would break) vs gates (propagate),
  * the union ancestor cone of the affected checks = the extended core,
  * the wide-multiply census of that cone -> QUBO size,
  * whether each affected check owns a private p-handle (self-absorbing).

Usage: python3 core_extend.py [state.json]
Writes: s13/core_extend.json
"""
import os, sys, json, time
from collections import deque, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

LAB = os.path.join(os.path.dirname(__file__), '..')
NA  = len(L.avars)
P   = 2**256 - 2**32 - 977
DRIVERS = [2081, 4287, 8731, 9118]


def degree(a):
    return max((len(m) for m in L.polys[a]), default=0)


def cone_of(atoms):
    need, seen = set(), set(atoms)
    dq = deque()
    for a in atoms:
        for x in L.avars[a]:
            if x not in need:
                need.add(x); dq.append(x)
    while dq:
        x = dq.popleft()
        d = L.definer.get(x)
        if d is None:
            continue
        seen.add(d)
        for y in L.avars[d]:
            if y not in need:
                need.add(y); dq.append(y)
    return need, seen


def main():
    t0 = time.time()
    state = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = L.load(state)
    free = {t for t in range(L.NVARS) if t not in L.definer}
    occ = defaultdict(list)
    for a in range(NA):
        for x in L.avars[a]:
            occ[x].append(a)
    solo = {x for x in free if len(occ[x]) == 1}

    print(f"[ext] state = {os.path.basename(state)}")
    print(f"[ext] drivers: {['x%d' % d for d in DRIVERS]}\n")

    # --- descendants: everything the drivers can disturb ---------------------
    users = defaultdict(list)                 # var -> atoms containing it
    for a in range(NA):
        for x in L.avars[a]:
            users[x].append(a)

    affected_checks, touched_atoms = set(), set()
    dq = deque(DRIVERS); seenv = set(DRIVERS)
    while dq:
        x = dq.popleft()
        for a in users[x]:
            touched_atoms.add(a)
            oc = L.atom_out.get(a)
            if oc is None:
                affected_checks.add(a)
            else:
                y = oc[1]
                if y not in seenv:
                    seenv.add(y); dq.append(y)
    print(f"[ext] driver descendants: {len(seenv):,} variables, "
          f"{len(touched_atoms):,} atoms touched")
    print(f"[ext] CHECKS at risk: {len(affected_checks):,}")

    # --- which of those checks self-absorb via a private p-handle ------------
    def near(a, depth=2):
        seen = set(L.avars[a]); fr = set(L.avars[a])
        for _ in range(depth):
            nx = set()
            for x in fr:
                d = L.definer.get(x)
                if d is not None:
                    nx |= L.avars[d]
            nx -= seen; seen |= nx; fr = nx
        return seen
    absorbing = {a for a in affected_checks if near(a) & solo}
    hard = sorted(affected_checks - absorbing)
    print(f"[ext]   of these, {len(absorbing):,} own a private p-handle "
          f"(self-absorbing)")
    print(f"[ext]   ** HARD checks (real constraints): {len(hard)} **")
    if hard[:24]:
        print(f"[ext]   hard: {['a%d' % a for a in hard[:24]]}"
              f"{' ...' if len(hard) > 24 else ''}")

    # --- extended core = cone of the hard checks + the 2 residual checks -----
    core_checks = sorted(set(hard) | {35759, 35760})
    cvars, catoms = cone_of(core_checks)
    cfree = sorted(cvars & free)
    print(f"\n[ext] EXTENDED CORE = cone of {len(core_checks)} checks")
    print(f"      {len(cvars):,} variables, {len(catoms):,} atoms, "
          f"{len(cfree):,} free inputs")

    # --- arithmetic census -> QUBO size --------------------------------------
    wide = {x for x in cvars if abs(v[x]).bit_length() >= 200}
    muls, small, lins = [], [], []
    for a in catoms:
        if degree(a) < 2:
            lins.append(a); continue
        exp = any(len(m) >= 2 and sum(1 for x in m if x in wide) >= 2
                  for m in L.polys[a])
        (muls if exp else small).append(a)
    print(f"      arithmetic: {len(lins):,} linear, {len(small):,} cheap products, "
          f"**{len(muls)} wide 256x256 multiplies**")

    LIMB = 16; nl = 256 // LIMB
    blk = LIMB * LIMB * 4                    # ~1,024 binary per limb block
    nblocks = len(muls) * nl * nl
    print(f"\n[ext] QUBO SIZE (16-bit limbs, {blk} binary per limb block)")
    print(f"      unknown residues : {len(cfree)} free x 256 bits "
          f"= {len(cfree)*256:,} binary")
    print(f"      multiply blocks  : {len(muls)} x {nl*nl} = {nblocks:,} blocks "
          f"-> {nblocks*blk:,} binary")
    print(f"      TOTAL            : ~{len(cfree)*256 + nblocks*blk:,} binary")
    print(f"      block size       : ~{blk} binary each "
          f"(target 1,000-5,000: {'OK' if 1000 <= blk <= 5000 else 'tune LIMB'})")
    print(f"      coupling         : carry words of ~{2*LIMB} bits between "
          f"adjacent limb blocks")

    json.dump({'drivers': DRIVERS, 'affected_checks': sorted(affected_checks),
               'absorbing': sorted(absorbing), 'hard_checks': hard,
               'core_checks': core_checks, 'core_vars': sorted(cvars),
               'core_atoms': sorted(catoms), 'core_free': cfree,
               'wide_muls': sorted(muls)},
              open(os.path.join(os.path.dirname(__file__), 'core_extend.json'), 'w'))
    print(f"\n[ext] {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
