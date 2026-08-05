#!/usr/bin/env python3
"""Is atom 40782 implied by the twist? Force x_18274:=x_9770, x_17728:=x_3183 in
the all-0 assignment and recompute 40782. Also classify each of 40782's 52 vars:
does it move with the 22 bits, the 233 bits, both, or neither?"""
import json, time
from confluent_eval5 import build5, make_forward
from propagate import atom_vars

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]

    base = solve(list(bestval), [])
    print(f"resid 40782 at all-0 = {resid(A[40782], base)}", flush=True)
    print(f"resid 1817 = {resid(A[1817], base)}", flush=True)
    print(f"resid 44271 = {resid(A[44271], base)}", flush=True)

    # force the twist
    v2 = list(base)
    v2[18274] = v2[9770]; v2[17728] = v2[3183]
    print(f"\nAfter forcing x_18274:=x_9770, x_17728:=x_3183:")
    print(f"  resid 1817 = {resid(A[1817], v2)}", flush=True)
    print(f"  resid 44271 = {resid(A[44271], v2)}", flush=True)
    print(f"  resid 30378 = {resid(A[30378], v2)}", flush=True)
    print(f"  resid 40782 = {resid(A[40782], v2)}", flush=True)

    # classify the 52 vars of atom 40782
    va = sorted(atom_vars(A[40782]))
    move22 = set(); move233 = set()
    # single-flip effect on each var
    b0 = {v: base[v] for v in va}
    for b in list(BITS22):
        val = solve(list(bestval), [b])
        for v in va:
            if val[v] != b0[v]: move22.add(v)
    for b in bits233[:40]:  # sample
        val = solve(list(bestval), [b])
        for v in va:
            if val[v] != b0[v]: move233.add(v)
    both = move22 & move233
    only22 = move22 - move233
    only233 = move233 - move22
    neither = set(va) - move22 - move233
    print(f"\natom 40782's 52 vars: only22={len(only22)} only233={len(only233)} both={len(both)} neither(const?)={len(neither)}", flush=True)
    print(f"  only22: {sorted(only22)}")
    print(f"  only233: {sorted(only233)}")
    print(f"  both: {sorted(both)}")
    print(f"  neither: {sorted(neither)}")
    # values of the 'neither' (constant) vars
    print("  neither values:", {v: base[v] for v in sorted(neither)})
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
