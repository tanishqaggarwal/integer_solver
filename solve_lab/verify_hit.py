#!/usr/bin/env python3
"""Verify a 22-bit hit pattern exactly in Z: forward-evaluate the whole circuit
with those bits set, then run the full checker."""
import json, sys, time
from collections import deque, defaultdict
from confluent_eval4 import build, boolean_vars
from propagate import load_atoms, NVARS

def forward_Z(setbits):
    atoms, gate_modP, loadinj, seq, _ = build()
    A = load_atoms()
    best = json.load(open('best/best_partial_39019.json'))
    val = [0]*NVARS
    for k, x in best.items(): val[int(k[2:])] = x
    prov = json.load(open('eval_order.json'))['prov']
    df = {}
    for v in range(NVARS):
        p = prov[v] if v < len(prov) else None
        if p and p[0] >= 0: df[v] = p[0]
    df[9770] = 27973; df[3183] = 27978
    bset = boolean_vars(A)
    load_atom = {}
    for a, poly in enumerate(A):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        big = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in bset and big and xB not in load_atom:
                load_atom[xB] = (a, bit); break
    for b in setbits: val[b] = 1
    ndiv = 0
    for v in seq:
        if v in loadinj:
            a, bit = load_atom[v]
            cbx = next(c for m, c in A[a].items() if len(m) == 2 and v in m and bit in m)
            if val[bit] == 0:
                val[v] = 0
            else:
                rest = 0
                for m, c in A[a].items():
                    if len(m) == 2 and v in m and bit in m: continue
                    t = c
                    for x in m: t *= (1 if x == bit else val[x])
                    rest += t
                num = -rest; den = cbx * 1
                if num % den != 0: ndiv += 1
                val[v] = num // den
        else:
            a = df[v]; coef = 0; rest = 0
            for m, c in A[a].items():
                if m == (v,): coef += c
                else:
                    t = c
                    for x in m: t *= val[x]
                    rest += t
            if coef == 0: continue
            if (-rest) % coef != 0: ndiv += 1
            val[v] = (-rest) // coef
    return val, ndiv, A

def main():
    hits = json.load(open('confluent_enum_hits.json'))
    print(f"{len(hits)} hit patterns to verify")
    for i, bits in enumerate(hits[:50]):
        val, ndiv, A = forward_Z(bits)
        viol = 0
        for poly in A:
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s: viol += 1
        print(f"  hit {i} bits={bits}: {viol} atoms violated, {ndiv} non-divisible", flush=True)
        if viol == 0:
            json.dump({f"x_{j}": val[j] for j in range(NVARS)}, open('cand_SOLVED.json', 'w'))
            print("  *** FULL SOLUTION! wrote cand_SOLVED.json ***", flush=True)
            return

if __name__ == '__main__':
    main()
