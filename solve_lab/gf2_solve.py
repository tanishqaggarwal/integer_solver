#!/usr/bin/env python3
"""Solve the system mod 2 by monomial-linearization + GF(2) Gaussian elimination.
Over GF(2) every var is a bit, x^2=x (multilinear), and a core bit's mod-2 value
equals its 0/1 value. If elimination forces core bits, we read them off directly."""
import json, time, sys
from collections import defaultdict

def monokey(m):
    # multilinear: distinct sorted vars
    return tuple(sorted(set(m)))

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    mc = json.load(open('main_comp.json'))
    mainv = set(mc['main_vars']); mainbits = mc['main_bits']

    # build GF(2) rows for atoms touching the main component
    mono_id = {}
    def mid(m):
        if m not in mono_id:
            mono_id[m] = len(mono_id)
        return mono_id[m]
    rows = []  # each row: (set of monomial-ids, const-bit)
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if not (vs & mainv):
            continue
        acc = defaultdict(int)
        for m, c in poly:
            acc[monokey(tuple(m))] ^= (c & 1)
        s = set(); const = 0
        for m, b in acc.items():
            if b:
                if m == ():
                    const ^= 1
                else:
                    s.add(mid(m))
        rows.append((s, const))
    print(f"GF(2): {len(rows)} rows, {len(mono_id)} monomials in {time.time()-t0:.1f}s", flush=True)

    # degree-1 monomial ids for core bits
    bit_mono = {}
    for b in mainbits:
        k = (b,)
        if k in mono_id:
            bit_mono[b] = mono_id[k]

    # sparse GF(2) Gaussian elimination
    pivot = {}   # monomial-id -> (row set, const)
    t1 = time.time()
    processed = 0
    for s, const in rows:
        s = set(s); c = const
        # reduce
        while True:
            red = False
            for mm in list(s):
                if mm in pivot:
                    ps, pc = pivot[mm]
                    s ^= ps; c ^= pc
                    red = True
                    break
            if not red:
                break
        if not s:
            continue  # 0 = c (consistent if c==0)
        p = min(s)
        pivot[p] = (frozenset(s), c)
        processed += 1
        if processed % 5000 == 0:
            print(f"   pivots={len(pivot)} ({time.time()-t1:.0f}s)", flush=True)
    print(f"rank {len(pivot)} / {len(mono_id)} monomials ({time.time()-t0:.0f}s)")

    # back-substitute to find fully-determined degree-1 vars (esp. core bits)
    # a monomial is 'determined' if its pivot row is just {itself} (=> equals const)
    forced = {}
    for p, (s, c) in pivot.items():
        if len(s) == 1:
            forced[p] = c
    # map forced monomials that are core bits
    forced_bits = {}
    for b, mm in bit_mono.items():
        # is bit b's monomial forced?
        if mm in forced:
            forced_bits[b] = forced[mm]
    print(f"core bits with a degree-1 monomial in system: {len(bit_mono)}")
    print(f"core bits FORCED mod 2 (singleton pivot): {len(forced_bits)}")
    if forced_bits:
        from collections import Counter
        print("  forced values:", dict(Counter(forced_bits.values())))
        print("  sample:", dict(list(forced_bits.items())[:20]))
    json.dump({str(k): v for k, v in forced_bits.items()}, open('gf2_forced_bits.json', 'w'))
    print(f"total time {time.time()-t0:.0f}s")

if __name__ == '__main__':
    main()
