#!/usr/bin/env python3
"""Identify structurally-relevant bit-pairs. A product x_b*x_c only becomes
nonzero when the bit loading x_b AND the bit loading x_c are both 1. So the
circuit is quadratic in bits with cross-terms only between such pairs. Compute
which wires each single bit activates, then collect the relevant pairs."""
import json, time
from collections import deque, defaultdict
import multiprocessing as mp
from propagate import load_atoms
from repair import boolean_vars
from modp import ModPEngine, P, NVARS

_G = {}
def init():
    atoms = load_atoms(); bset = boolean_vars(atoms)
    base = ModPEngine(atoms); base.propagate()
    _G.update(atoms=atoms, bset=bset, bval=list(base.val))

def activated_wires(bit):
    atoms = _G['atoms']; bset = _G['bset']
    eng = ModPEngine(atoms); eng.val = list(_G['bval']); eng.wl = deque(); eng.inwl = [False]*len(atoms)
    if eng.val[bit] is None: eng.assign(bit, 1)
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    # wires nonzero (that were zero at baseline)
    nz = [v for v in range(NVARS) if eng.val[v] not in (None, 0) and _G['bval'][v] in (None, 0)]
    return bit, nz

def main():
    t0 = time.time()
    control = json.load(open('control_bits.json'))
    atoms = load_atoms(); mainv = set(json.load(open('main_comp.json'))['main_vars'])
    # wire -> set of bits that activate it
    wire_bits = defaultdict(set)
    with mp.Pool(4, initializer=init) as pool:
        for k, (bit, nz) in enumerate(pool.imap_unordered(activated_wires, control)):
            for w in nz: wire_bits[w].add(bit)
            if (k+1) % 64 == 0: print(f"  {k+1}/{len(control)} bits, {time.time()-t0:.0f}s", flush=True)
    print(f"wires activated by >=1 bit: {len(wire_bits)} ({time.time()-t0:.0f}s)", flush=True)
    hist = defaultdict(int)
    for w, bs in wire_bits.items(): hist[len(bs)] += 1
    print("  #bits-activating-a-wire histogram:", dict(sorted(hist.items())[:12]))

    # relevant pairs from product monomials
    pairs = defaultdict(int)   # (b1,b2) -> count of products it enables
    singleton_products = 0
    for poly in atoms:
        vs = set()
        for m in poly: vs.update(m)
        if not (vs & mainv): continue
        for m in poly:
            if len(m) >= 2:
                # operands and their activating bits
                bitsets = [wire_bits.get(u, set()) for u in m]
                # relevant cross pairs
                allbits = set()
                for bs in bitsets: allbits |= bs
                # pairs across different operands
                for i in range(len(m)):
                    for j in range(i+1, len(m)):
                        for a in bitsets[i]:
                            for b in bitsets[j]:
                                if a != b: pairs[tuple(sorted((a, b)))] += 1
    print(f"distinct relevant bit-pairs: {len(pairs)} ({time.time()-t0:.0f}s)")
    top = sorted(pairs.items(), key=lambda kv: -kv[1])[:10]
    print("  top pairs by #products enabled:", top)
    json.dump({f"{a}_{b}": c for (a, b), c in pairs.items()}, open('relevant_pairs.json', 'w'))
    json.dump({str(w): sorted(bs) for w, bs in wire_bits.items()}, open('wire_bits.json', 'w'))
    print("wrote relevant_pairs.json, wire_bits.json")

if __name__ == '__main__':
    main()
