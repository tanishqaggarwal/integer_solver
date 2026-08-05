#!/usr/bin/env python3
"""Confluent forward evaluator WITH residue loading (the missing mechanism).

Huge-atoms bit*(x_B - HUGE) - s*x_C = 0 make x_B a bit-gated value-input:
x_B = 0 if bit=0, else HUGE + s*x_C. In best prov x_B is defined by a *different*
gate (producing 0, since bit=0), so a naive evaluator never loads residues when a
bit is flipped. Here we override every huge-atom target to load conditionally on
its bit, giving a correct confluent forward evaluation. Then screen bits (mod-P):
now a bit flip actually loads a residue and propagates."""
import json, time, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv
sys.setrecursionlimit(1_000_000)

def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        if len(poly) == 2:
            ms = sorted(poly.keys(), key=len)
            if len(ms[0]) == 1 and len(ms[1]) == 2 and ms[1] == (ms[0][0], ms[0][0]) and poly[ms[0]] == -poly[ms[1]]:
                bset.add(ms[0][0])
    return bset

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    bset = boolean_vars(atoms)
    # detect huge-atoms: a monomial (bit, x_B) with bit boolean and a huge coef somewhere
    load = {}   # x_B -> (atom_idx, bit, coef_bxB)
    for a, poly in enumerate(atoms):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        # which of the two is the bit?
        cand = [(m2[0], m2[1]), (m2[1], m2[0])]
        bigcoef = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in cand:
            if bit in bset and bigcoef and xB not in load:
                load[xB] = (a, bit, c2 % P)
                break
    print(f"huge-atom loads detected: {len(load)} ({time.time()-t0:.0f}s)", flush=True)

    prov = json.load(open('eval_order.json'))['prov']
    order = json.load(open('eval_order.json'))['order']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x

    df = [None]*NVARS
    for v, p in enumerate(prov):
        if p and p[0] >= 0: df[v] = p[0]
    df[9770] = 27973; df[3183] = 27978

    # precompile normal gates; load vars handled specially
    gate = {}
    for v in range(NVARS):
        if v in load: continue
        a = df[v]
        if a is None: continue
        coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef = (coef + c) % P
            else: terms.append((c % P, m))
        gate[v] = (coef, terms)
    # load-var terms: rest monomials (all except the (bit,x_B) one)
    loadinfo = {}
    for xB, (a, bit, cbx) in load.items():
        terms = [(c % P, m) for m, c in atoms[a].items() if not (len(m) == 2 and xB in m and bit in m)]
        loadinfo[xB] = (bit, cbx, terms)

    bm = [x % P for x in bestval]
    # sequence: eval order; ensure load vars & derived appear
    seqset = set(gate) | set(load)
    seq = [v for v in order if v in seqset]
    seq += [v for v in seqset if v not in set(order)]

    CHECKS = [1817, 30378, 44271]
    cp = [[(c % P, m) for m, c in atoms[a].items()] for a in CHECKS]

    def forward(setbits):
        val = list(bm)
        sb = set(setbits)
        for b in setbits: val[b] = 1
        for v in seq:
            if v in loadinfo:
                bit, cbx, terms = loadinfo[v]
                if val[bit] == 0:
                    val[v] = 0
                else:
                    rs = 0
                    for c, m in terms:
                        t = c
                        for x in m: t = (t * val[x]) % P
                        rs = (rs + t) % P
                    val[v] = (-rs * inv((cbx * val[bit]) % P)) % P
            else:
                coef, terms = gate[v]
                if coef == 0: continue
                rs = 0
                for c, m in terms:
                    t = c
                    for x in m: t = (t * val[x]) % P
                    rs = (rs + t) % P
                val[v] = (-rs * inv(coef)) % P
        res = []
        for poly in cp:
            s = 0
            for c, m in poly:
                t = c
                for x in m: t = (t * val[x]) % P
                s = (s + t) % P
            res.append(s)
        return res, val

    base, _ = forward([])
    print(f"baseline checks (0=pass): {[1 if x else 0 for x in base]} ({time.time()-t0:.0f}s)", flush=True)

    control = json.load(open('control_bits.json'))
    changed = []
    for b in control:
        r, _ = forward([b])
        if r != base:
            npass = sum(1 for x in r if x == 0)
            changed.append((npass, b, [1 if x else 0 for x in r]))
    changed.sort(reverse=True)
    print(f"single bits that CHANGE the checks: {len(changed)} ({time.time()-t0:.0f}s)", flush=True)
    for npass, b, r in changed[:25]:
        print(f"  bit {b}: pass {npass}/3  {r}", flush=True)
    json.dump([[n, b, r] for n, b, r in changed], open('confluent2_single.json', 'w'))
    print("wrote confluent2_single.json", flush=True)

if __name__ == '__main__':
    main()
