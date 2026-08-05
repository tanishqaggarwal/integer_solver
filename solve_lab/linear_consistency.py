#!/usr/bin/env python3
"""Exact GF(p) Gaussian elimination on all LINEAR atoms to test whether the
'zero every atom' system's linear part is even consistent, and to find its rank
and which variables (esp. the 256 core bits) it forces."""
import json, time, sys
from collections import defaultdict

P = (1 << 61) - 1  # Mersenne prime

def inv(a): return pow(a % P, P - 2, P)

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    # linear atoms: degree exactly 1
    rows = []
    for poly in atoms:
        deg = max((len(m) for m, c in poly), default=0)
        if deg != 1:
            continue
        row = {}
        for m, c in poly:
            key = m[0] if len(m) == 1 else 'c'
            row[key] = (row.get(key, 0) + c) % P
        rows.append(row)
    print(f"linear atoms: {len(rows)}  ({time.time()-t0:.1f}s)")

    # sparse Gaussian elimination
    pivot = {}   # var -> reduced row (dict) with that var as pivot (coef 1)
    inconsistent = 0
    order = sorted(range(len(rows)), key=lambda i: len(rows[i]))
    t1 = time.time()
    for cnt, i in enumerate(order):
        row = dict(rows[i])
        # reduce by existing pivots
        changed = True
        while changed:
            changed = False
            for v in [k for k in row if k != 'c']:
                if v in pivot and row.get(v, 0) != 0:
                    f = row[v]
                    pr = pivot[v]
                    for k, cf in pr.items():
                        row[k] = (row.get(k, 0) - f * cf) % P
                    if row.get(v):
                        row[v] = 0
                    row = {k: c for k, c in row.items() if c % P != 0}
                    changed = True
                    break
        vars_left = [k for k in row if k != 'c']
        if not vars_left:
            if row.get('c', 0) % P != 0:
                inconsistent += 1
            continue
        # new pivot on the smallest-index var
        pv = min(vars_left)
        f = inv(row[pv])
        prow = {k: (c * f) % P for k, c in row.items()}
        prow[pv] = 1
        pivot[pv] = prow
        if (cnt + 1) % 4000 == 0:
            print(f"   processed {cnt+1}/{len(rows)}, pivots={len(pivot)}, incons={inconsistent}, {time.time()-t1:.0f}s", flush=True)

    print(f"\nrank (pivots): {len(pivot)}")
    print(f"inconsistent rows (0 = nonzero): {inconsistent}")
    # which core bits are forced to a constant (pivot row with only 'c')?
    mainbits = set(json.load(open('main_comp.json'))['main_bits'])
    forced_bits = {}
    for b in mainbits:
        if b in pivot:
            pr = pivot[b]
            others = [k for k in pr if k not in ('c', b)]
            if not others:
                forced_bits[b] = (-pr.get('c', 0) * inv(pr[b])) % P
    print(f"core bits appearing as pivots: {sum(1 for b in mainbits if b in pivot)}")
    print(f"core bits FORCED to a constant by linear system: {len(forced_bits)}")
    if forced_bits:
        print("  sample forced bits:", dict(list(forced_bits.items())[:20]))
    print(f"total time {time.time()-t0:.0f}s")
    json.dump({"rank": len(pivot), "inconsistent": inconsistent,
               "forced_core_bits": {str(k): int(v) for k, v in forced_bits.items()}},
              open('linear_consistency.json', 'w'))

if __name__ == '__main__':
    main()
