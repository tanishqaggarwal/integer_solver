#!/usr/bin/env python3
"""Custom decompose-and-solve: atoms 27973/27978/45004 depend only on the same
22 bits. Build a fast mod-P evaluator for the sub-circuit those 22 bits drive,
enumerate all 2^22 patterns, and collect the ones that zero all three atoms.
Then check each survivor against atom 41470 and the full system over Z."""
import json, time, itertools
from collections import deque
from propagate import load_atoms
from repair import boolean_vars
from modp import ModPEngine, P, NVARS, inv

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
TARGETS = [27973, 27978, 45004]

def main():
    t0 = time.time()
    atoms = load_atoms(); bset = boolean_vars(atoms)
    # baseline (all bits 0) mod P
    base = ModPEngine(atoms); base.propagate()
    for b in [x for x in bset if base.val[x] is None]:
        if base.val[b] is None: base.assign(b, 0); base.propagate()
    for v in range(NVARS):
        if base.val[v] is None: base.assign(v, 0); base.propagate()
    base_val = [x if x is not None else 0 for x in base.val]

    # cone of the 22 bits: set them all =1, see which wires change
    eng = ModPEngine(atoms); eng.val = list(base.val)  # careful: reuse determined
    # rebuild from scratch with 22 bits = 1
    e2 = ModPEngine(atoms)
    for b in BITS22: e2.assign(b, 1)
    e2.propagate()
    for b in [x for x in bset if e2.val[x] is None]:
        if e2.val[b] is None: e2.assign(b, 0); e2.propagate()
    for v in range(NVARS):
        if e2.val[v] is None: e2.assign(v, 0); e2.propagate()
    allon = [x if x is not None else 0 for x in e2.val]
    cone = [v for v in range(NVARS) if allon[v] != base_val[v]]
    print(f"cone of 22 bits: {len(cone)} wires change ({time.time()-t0:.0f}s)", flush=True)

    # We need each cone wire's defining atom to recompute it. Use provenance order.
    eo = json.load(open('eval_order.json')); order = eo['order']; prov = eo['prov']
    coneset = set(cone) | set(BITS22)
    # sub-order: cone wires in eval order, with their defining atoms
    sub = [(v, prov[v]) for v in order if v in coneset and prov[v] and prov[v][0] >= 0]
    defatom = {v: [(c % P, tuple(m)) for m, c in atoms[p[0]].items()] for v, p in sub}
    suborder = [v for v, p in sub]
    print(f"sub-circuit recompute wires: {len(suborder)}", flush=True)
    # target atoms as monomial lists
    tgt = {ai: [(c % P, tuple(m)) for m, c in atoms[ai].items()] for ai in TARGETS}

    val = list(base_val)
    def evalpat(pat):
        # pat: tuple of 0/1 for BITS22
        for i, b in enumerate(BITS22): val[b] = pat[i]
        for v in suborder:
            c0 = c1 = 0
            for coef, m in defatom[v]:
                if v in m:
                    t = coef
                    for x in m:
                        if x != v: t = (t*val[x]) % P
                    c1 = (c1 + t) % P
                else:
                    t = coef
                    for x in m: t = (t*val[x]) % P
                    c0 = (c0 + t) % P
            val[v] = (-c0*inv(c1)) % P if c1 else 0
        res = []
        for ai in TARGETS:
            s = 0
            for coef, m in tgt[ai]:
                t = coef
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            res.append(s)
        return res

    # sanity: baseline pattern
    r0 = evalpat(tuple([0]*22))
    print(f"baseline (22 bits=0) target atom values mod P: {[1 if x else 0 for x in r0]} (nonzero=violated)", flush=True)
    t1 = time.time()
    hits = []
    N = 1 << 22
    for code in range(N):
        pat = tuple((code >> i) & 1 for i in range(22))
        r = evalpat(pat)
        if not any(r):
            hits.append(pat)
            print(f"  HIT pattern {code:#08x}: {[b for i,b in enumerate(BITS22) if pat[i]]}", flush=True)
        if code % 500000 == 0 and code:
            print(f"  ...{code}/{N} ({time.time()-t1:.0f}s), hits={len(hits)}", flush=True)
    print(f"enumeration done {time.time()-t1:.0f}s; {len(hits)} patterns zero all 3 atoms", flush=True)
    json.dump([[b for i,b in enumerate(BITS22) if pat[i]] for pat in hits], open('sub_hits.json','w'))

if __name__ == '__main__':
    main()
