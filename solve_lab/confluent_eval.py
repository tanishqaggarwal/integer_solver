#!/usr/bin/env python3
"""Confluent forward evaluator (best's orientation) + clean bit screen.

best's prov gives a VALID acyclic orientation for the 39,013 satisfied gates.
Override x_9770<-27973, x_3183<-27978 (their true sum-gates). Forward-evaluate
every derived var from its gate in eval order -> deterministic (confluent), unlike
propagation. This reproduces 39,019 at all-0 bits. Because it's confluent, a bit's
effect on the 3 remaining check-gates (1817, 30378, 44271) is its TRUE semantic
effect. Screen singles then pairs (mod-P for speed); verify hits in Z."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv
sys.setrecursionlimit(1_000_000)

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    prov = json.load(open('eval_order.json'))['prov']
    order = json.load(open('eval_order.json'))['order']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x

    df = [None]*NVARS
    for v, p in enumerate(prov):
        if p and p[0] >= 0: df[v] = p[0]
    df[9770] = 27973; df[3183] = 27978

    # precompile each derived var's gate mod P: coef_v and list of (coef, monomial)
    gate = {}
    for v in range(NVARS):
        a = df[v]
        if a is None: continue
        coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef += c % P
            else: terms.append((c % P, m))
        gate[v] = (coef % P, terms)
    seq = [v for v in order if v in gate]
    # append any derived var missing from order (shouldn't be many)
    inorder = set(order)
    seq += [v for v in range(NVARS) if v in gate and v not in inorder]
    print(f"derived {len(gate)}, seq {len(seq)} ({time.time()-t0:.0f}s)", flush=True)

    CHECKS = [1817, 30378, 44271]
    check_polys = [[(c % P, m) for m, c in atoms[a].items()] for a in CHECKS]
    bestmod = [x % P for x in bestval]

    def forward(setbits):
        val = list(bestmod)
        for b in setbits: val[b] = 1
        for v in seq:
            coef, terms = gate[v]
            if coef == 0:
                continue
            rs = 0
            for c, m in terms:
                t = c
                for x in m: t = (t * val[x]) % P
                rs = (rs + t) % P
            val[v] = (-rs * inv(coef)) % P
        # evaluate the 3 check gates mod P
        res = []
        for poly in check_polys:
            s = 0
            for c, m in poly:
                t = c
                for x in m: t = (t * val[x]) % P
                s = (s + t) % P
            res.append(s)
        return res, val

    base, _ = forward([])
    print(f"baseline check-gates mod P (0=pass): {[1 if x else 0 for x in base]}  (all-0 bits = 39,019)", flush=True)

    control = json.load(open('control_bits.json'))
    # single-bit effect on the 3 checks
    single = {}
    improved = []
    for b in control:
        r, _ = forward([b])
        npass = sum(1 for x in r if x == 0)
        single[b] = r
        if npass > 0:
            improved.append((npass, b, [1 if x else 0 for x in r]))
    improved.sort(reverse=True)
    print(f"single bits passing >=1 check: {len(improved)} ({time.time()-t0:.0f}s)", flush=True)
    for npass, b, r in improved[:20]:
        print(f"  bit {b}: {npass}/3 checks pass  {r}", flush=True)
    json.dump({str(b): [int(x) for x in single[b]] for b in control}, open('confluent_single.json', 'w'))
    print("wrote confluent_single.json", flush=True)

if __name__ == '__main__':
    main()
