#!/usr/bin/env python3
"""
Measure the ARITHMETIC of the core: the ancestor cone of the binding congruences,
and how many 256-bit modular multiplies a QUBO of it must contain.

Fixes the R5 test of core_reduce.py: a p-handle is rarely inside the check atom
itself -- it sits one gate up (the check references a gate variable whose definer
is `wire*h`).  Here a check counts as SELF-ABSORBING if any variable in its cone,
at depth <= 2, is a free input occurring in exactly one atom.

Then, for the checks that remain binding at the given state (the true residual),
we walk the ancestor cone and count:
    * cone size (variables / atoms),
    * MUL = atoms of degree >= 2 whose operands are both 256-bit-wide values
      (these are the expensive modular multiplies),
    * LIN = atoms that are linear (free in a QUBO: they fold into the objective).

QUBO cost model at 256-bit width, 16-bit limbs:
    one 256x256 schoolbook multiply = 16*16 = 256 limb-products,
    each limb-product = a 16x16 multiply block ~ 16*16 partial-product bits
      + accumulation ~= 1,024 binary variables  -> a 1k-block,
    so one modular multiply ~= 256 blocks ~= 262,144 binary variables.
Limb blocks couple ONLY through carry/accumulator words (~32 bits), which is the
minimal-coupling structure this decomposition targets.

Usage: python3 core_cone.py [state.json]
"""
import os, sys, json, time
from collections import defaultdict, deque
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "s10"))
import ad

LAB = os.path.join(os.path.dirname(__file__), '..')
NA  = len(L.avars)
P   = 2**256 - 2**32 - 977


def degree(a):
    return max((len(m) for m in L.polys[a]), default=0)


def main():
    t0 = time.time()
    state = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = L.load(state)
    print(f"[cone] state = {os.path.basename(state)}")

    free = {t for t in range(L.NVARS) if t not in L.definer}
    occ = defaultdict(int)
    for a in range(NA):
        for x in L.avars[a]:
            occ[x] += 1
    solo = {x for x in free if occ[x] == 1}

    checks = [a for a in range(NA) if L.atom_out.get(a) is None]

    # --- corrected self-absorption test: solo handle within depth 2 ----------
    def near_vars(a, depth=2):
        seen = set(L.avars[a]); frontier = set(L.avars[a])
        for _ in range(depth):
            nxt = set()
            for x in frontier:
                d = L.definer.get(x)
                if d is not None:
                    nxt |= L.avars[d]
            nxt -= seen; seen |= nxt; frontier = nxt
        return seen

    absorbed = [a for a in checks if near_vars(a) & solo]
    binding  = [a for a in checks if a not in set(absorbed)]
    print(f"[cone] R5 (depth-2 handle test): {len(absorbed):,} checks self-absorb, "
          f"{len(binding):,} binding")

    # --- the TRUE residual at this state: checks that are nonzero ------------
    av = L.all_atom_values(v)
    nz = [a for a in checks if av[a] != 0]
    print(f"[cone] nonzero checks at this state: {len(nz)} -> {nz}")

    # --- ancestor cone of the nonzero checks ---------------------------------
    need, atoms_in = set(), set()
    dq = deque()
    for a in nz:
        atoms_in.add(a)
        for x in L.avars[a]:
            if x not in need:
                need.add(x); dq.append(x)
    while dq:
        x = dq.popleft()
        d = L.definer.get(x)
        if d is None:
            continue
        atoms_in.add(d)
        for y in L.avars[d]:
            if y not in need:
                need.add(y); dq.append(y)
    cone_free = sorted(need & free)
    print(f"[cone] ancestor cone: {len(need):,} variables, {len(atoms_in):,} atoms, "
          f"{len(cone_free):,} free inputs")

    # --- arithmetic census over the cone -------------------------------------
    wide = {x for x in need if v[x].bit_length() >= 200}
    muls, lins, small = [], [], []
    for a in atoms_in:
        d = degree(a)
        if d < 2:
            lins.append(a); continue
        # a multiply is expensive iff BOTH operands are wide
        expensive = False
        for mono in L.polys[a]:
            if len(mono) >= 2 and sum(1 for x in mono if x in wide) >= 2:
                expensive = True; break
        (muls if expensive else small).append(a)
    print(f"[cone] arithmetic: {len(lins):,} linear atoms (free in QUBO), "
          f"{len(small):,} cheap products (one wide operand), "
          f"**{len(muls)} WIDE 256x256 multiplies**")
    if muls:
        print(f"[cone] wide-multiply atoms: {muls[:20]}")

    # --- QUBO sizing ----------------------------------------------------------
    LIMB = 16
    nlimb = 256 // LIMB
    per_limb_block = LIMB * LIMB * 4      # partial products + accumulate ~1,024
    per_mul_blocks = nlimb * nlimb        # 256 limb-products
    print(f"\n[cone] QUBO SIZING (16-bit limbs)")
    print(f"  one 256x256 modular multiply = {per_mul_blocks} limb blocks "
          f"x ~{per_limb_block} binary = {per_mul_blocks*per_limb_block:,} binary")
    print(f"  cone needs {len(muls)} wide multiplies "
          f"-> {len(muls)*per_mul_blocks:,} blocks, "
          f"{len(muls)*per_mul_blocks*per_limb_block:,} binary total")
    print(f"  advice unknowns: 13 x 256 = 3,328 binary")
    print(f"  => blocks of ~{per_limb_block} binary each, coupled only by "
          f"carry words (~{2*LIMB} bits per boundary)")

    json.dump({'binding': binding, 'nonzero_checks': nz,
               'cone_vars': sorted(need), 'cone_atoms': sorted(atoms_in),
               'cone_free': cone_free, 'wide_muls': muls},
              open(os.path.join(os.path.dirname(__file__), 'core_cone.json'), 'w'))
    print(f"\n[cone] {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
