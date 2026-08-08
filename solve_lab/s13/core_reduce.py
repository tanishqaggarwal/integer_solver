#!/usr/bin/env python3
"""
Reduce the instance to its CORE, then size the core for QUBO decomposition.

Chain of reductions, each measured (not assumed):

  R0  raw file                     39,033 equations / 38,748 variables
  R1  shared sub-expressions       42,267 atoms
  R2  circuit orientation (s9)     31,475 gates (DEFINED vars, functions of the
                                   free inputs) + 7,273 free inputs + checks.
                                   A gate is not an unknown: it is evaluated.
  R3  live free inputs             free inputs that actually reach a CHECK
  R4  the advice                   free inputs carrying real information
                                   (~296-bit numbers) + the real message bits
  R5  the binding checks           checks whose value is NOT absorbed by a
                                   private p-handle -- see below.

The R5 test is the important one.  Nearly every check in this instance has the
shape   c*(A - B) - wire*h = 0   with `wire` a variable equal to p and `h` a free
input occurring in NO other atom (a "solo handle").  Over Z such a check is
satisfiable for ANY A,B by choosing h = c(A-B)/p -- PROVIDED that quotient is an
integer.  So the check contributes no equation to the core; it contributes the
divisibility  p | c*(A-B).  Checks with no private handle are hard constraints.

Usage:  python3 core_reduce.py [state.json]
Writes: s13/core.json  (core variable / check manifest)
"""
import os, sys, json, time
from collections import defaultdict, deque
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

LAB  = os.path.join(os.path.dirname(__file__), '..')
OUT  = os.path.join(os.path.dirname(__file__), 'core.json')
NA   = len(L.avars)
P    = 2**256 - 2**32 - 977


def main():
    t0 = time.time()
    state = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = L.load(state)
    print(f"[core] state = {os.path.basename(state)}")

    free  = [t for t in range(L.NVARS) if t not in L.definer]
    freeset = set(free)
    gates = [a for a in range(NA) if L.atom_out.get(a) is not None]
    checks = [a for a in range(NA) if L.atom_out.get(a) is None]
    print(f"[core] R2: {len(gates):,} gates, {len(free):,} free inputs, "
          f"{len(checks):,} checks")

    # ---- R5 first: which checks own a PRIVATE solo handle -------------------
    occ = defaultdict(int)                       # var -> #atoms containing it
    for a in range(NA):
        for x in L.avars[a]:
            occ[x] += 1
    solo_handles = {x for x in free if occ[x] == 1}
    absorbed, binding = [], []
    for a in checks:
        if any((x in solo_handles) for x in L.avars[a]):
            absorbed.append(a)
        else:
            binding.append(a)
    print(f"[core] R5: {len(solo_handles):,} solo free handles -> "
          f"{len(absorbed):,} checks self-absorb, "
          f"**{len(binding):,} checks are BINDING**")

    # ---- R4: classify free inputs by the information they carry -------------
    bits, advice, zero, other = [], [], [], []
    for x in free:
        b = v[x].bit_length()
        if v[x] in (0, 1) and b <= 1:
            (bits if v[x] == 1 else zero).append(x)
        elif 280 <= b <= 330:
            advice.append(x)
        else:
            other.append(x)
    print(f"[core] R4: free inputs -> {len(advice)} advice (~296-bit), "
          f"{len(bits):,} set-bits, {len(zero):,} zero, {len(other):,} other")
    print(f"[core]     advice vars: {[f'x{t}' for t in advice]}")

    # ---- R3: which free inputs actually reach a BINDING check ---------------
    # forward reachability in the DAG from each free input is expensive; instead
    # walk backwards from the binding checks over definer edges.
    need = set()
    dq = deque()
    for a in binding:
        for x in L.avars[a]:
            if x not in need:
                need.add(x); dq.append(x)
    while dq:
        x = dq.popleft()
        d = L.definer.get(x)
        if d is None:
            continue
        for y in L.avars[d]:
            if y not in need:
                need.add(y); dq.append(y)
    live_free = sorted(freeset & need)
    print(f"[core] R3: cone of the binding checks = {len(need):,} variables, "
          f"of which **{len(live_free):,} are free inputs**")

    live_advice = [x for x in live_free if x in set(advice)]
    live_bits   = [x for x in live_free if x in set(bits)]
    print(f"[core]     live advice = {len(live_advice)}  "
          f"{[f'x{t}' for t in live_advice]}")
    print(f"[core]     live set-bits = {len(live_bits)}")

    # ---- the core unknown budget -------------------------------------------
    # advice values are 296-bit; only their residue mod p matters (256 bits)
    # plus ~40 bits of k*p slack that the handles absorb.
    n_adv = len(live_advice)
    print(f"\n[core] CORE UNKNOWN BUDGET")
    print(f"  advice residues   : {n_adv} x 256 bits = {n_adv*256:,} binary")
    print(f"  message bits      : {len(live_bits):,} binary")
    print(f"  TOTAL genuine unknowns ~= {n_adv*256 + len(live_bits):,} binary")

    # ---- binding-check arithmetic weight ------------------------------------
    def degree(a):
        return max((len(m) for m in L.polys[a]), default=0)
    dh = defaultdict(int)
    for a in binding:
        dh[degree(a)] += 1
    print(f"  binding checks by degree: {dict(sorted(dh.items()))}")

    json.dump({
        'state': os.path.basename(state),
        'n_gates': len(gates), 'n_free': len(free), 'n_checks': len(checks),
        'solo_handles': len(solo_handles),
        'absorbed_checks': len(absorbed), 'binding_checks': binding,
        'advice': advice, 'live_advice': live_advice,
        'live_free': live_free, 'live_bits': live_bits,
        'cone_vars': sorted(need),
    }, open(OUT, 'w'))
    print(f"\n[core] manifest -> {OUT}   ({time.time()-t0:.1f}s)")


if __name__ == '__main__':
    main()
