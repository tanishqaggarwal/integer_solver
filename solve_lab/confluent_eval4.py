#!/usr/bin/env python3
"""Confluent evaluator via VALUE INJECTION on best's acyclic order (no cycles).

Use best's eval_order (acyclic) as-is; move x_9770/x_3183 to the end so their
sum-gate inputs are ready. For each huge-atom load-target x_B: after its best-prov
gate computes it, if the bit is 1 inject the loaded residue x_B = HUGE + s*x_C.
This models residue selection by free bits without re-orienting (hence no cycles),
giving a clean deterministic evaluation. Screen bits mod-P; verify hits in Z."""
import json, time, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv
sys.setrecursionlimit(2_000_000)

def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        if len(poly) == 2:
            ms = sorted(poly.keys(), key=len)
            if len(ms[0]) == 1 and len(ms[1]) == 2 and ms[1] == (ms[0][0], ms[0][0]) and poly[ms[0]] == -poly[ms[1]]:
                bset.add(ms[0][0])
    return bset

def build():
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    bset = boolean_vars(atoms)
    load = {}
    for a, poly in enumerate(atoms):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        bigcoef = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in bset and bigcoef and xB not in load:
                load[xB] = (a, bit, c2 % P); break
    prov = json.load(open('eval_order.json'))['prov']
    order = json.load(open('eval_order.json'))['order']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x

    df = [None]*NVARS
    for v, p in enumerate(prov):
        if p and p[0] >= 0: df[v] = p[0]
    # gates
    gate = {}
    for v in range(NVARS):
        a = df[v]
        if a is None: continue
        coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef = (coef + c) % P
            else: terms.append((c % P, m))
        gate[v] = (coef, terms)
    # override gate for 9770/3183 to their sum-gates
    for v, a in [(9770, 27973), (3183, 27978)]:
        coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef = (coef + c) % P
            else: terms.append((c % P, m))
        gate[v] = (coef, terms)
    # load injection terms (all monomials except the bit*xB product)
    loadinj = {}
    for xB, (a, bit, cbx) in load.items():
        terms = [(c % P, m) for m, c in atoms[a].items() if not (len(m) == 2 and xB in m and bit in m)]
        loadinj[xB] = (bit, cbx, terms)

    # order: best eval_order, but move 9770/3183 to the very end
    seq = [v for v in order if v in gate and v not in (9770, 3183)]
    tail = [v for v in (9770, 3183) if v in gate]
    seq = seq + tail
    # append any gate var missing from order
    inseq = set(seq)
    seq += [v for v in gate if v not in inseq]
    return atoms, gate, loadinj, seq, [x % P for x in bestval]

def main():
    t0 = time.time()
    atoms, gate, loadinj, seq, bm = build()
    print(f"built: {len(seq)} derived, {len(loadinj)} loads ({time.time()-t0:.0f}s)", flush=True)
    CHECKS = [1817, 30378, 44271]
    cp = [[(c % P, m) for m, c in atoms[a].items()] for a in CHECKS]

    def forward(setbits):
        val = list(bm)
        sb = set(setbits)
        for b in setbits: val[b] = 1
        for v in seq:
            coef, terms = gate[v]
            if coef != 0:
                rs = 0
                for c, m in terms:
                    t = c
                    for x in m: t = (t * val[x]) % P
                    rs = (rs + t) % P
                val[v] = (-rs * inv(coef)) % P
            # load injection
            inj = loadinj.get(v)
            if inj is not None:
                bit, cbx, lterms = inj
                if val[bit] != 0:
                    rs = 0
                    for c, m in lterms:
                        t = c
                        for x in m: t = (t * val[x]) % P
                        rs = (rs + t) % P
                    val[v] = (-rs * inv((cbx * val[bit]) % P)) % P
                else:
                    val[v] = 0
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
    print(f"single bits that CHANGE checks: {len(changed)} ({time.time()-t0:.0f}s)", flush=True)
    for npass, b, r in changed[:30]:
        print(f"  bit {b}: pass {npass}/3  {r}", flush=True)
    json.dump([[n, b, r] for n, b, r in changed], open('confluent4_single.json', 'w'))

if __name__ == '__main__':
    main()
