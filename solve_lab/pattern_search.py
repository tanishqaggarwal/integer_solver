#!/usr/bin/env python3
"""Vectorized search over arbitrary bit-patterns (confluent evaluator).

Evaluates a supplied list of bit-subsets against the two twist checks
(x_18274=x_9770, x_17728=x_3183) using the validated confluent evaluator over the
checks' backward cone, mod two 31-bit primes. Collects patterns zeroing both,
then verifies exactly in Z. Used to search pairs / triples that activate the
product terms the 22 linear-effect bits cannot reach alone."""
import json, time, sys, itertools
import numpy as np
from collections import deque, defaultdict
from confluent_eval4 import build, boolean_vars
from propagate import load_atoms, atom_vars, NVARS

CHECKS = [1817, 30378, 44271]
PRIMES = [2147483647, 2147483629]
BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def build_structures():
    atoms, gate, loadinj, seq, _ = build()
    A = load_atoms()
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x
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
    # backward cone of check vars
    deps = defaultdict(set)
    for v in gate:
        for c, m in gate[v][1]: deps[v].update(m)
        if v in loadinj:
            bit, cbx, lt = loadinj[v]; deps[v].add(bit)
            for c, m in lt: deps[v].update(m)
        deps[v].discard(v)
    checkvars = set()
    for a in CHECKS: checkvars.update(atom_vars(A[a]))
    seen = set(checkvars); dq = deque(checkvars)
    while dq:
        v = dq.popleft()
        for x in deps.get(v, ()):
            if x not in seen: seen.add(x); dq.append(x)
    conewires = [v for v in seq if v in seen and v in gate]
    exgate = {}; exload = {}
    for v in conewires:
        if v in loadinj:
            a, bit = load_atom[v]
            cbx = next(c for m, c in A[a].items() if len(m) == 2 and v in m and bit in m)
            lterms = [(c, m) for m, c in A[a].items() if not (len(m) == 2 and v in m and bit in m)]
            exload[v] = (bit, cbx, lterms)
        else:
            a = df[v]; coef = 0; terms = []
            for m, c in A[a].items():
                if m == (v,): coef += c
                else: terms.append((c, m))
            exgate[v] = (coef, terms)
    checkpoly = {a: list(A[a].items()) for a in CHECKS}
    return conewires, exgate, exload, checkpoly, bestval

def search(patterns, tag):
    t0 = time.time()
    conewires, exgate, exload, checkpoly, bestval = build_structures()
    print(f"[{tag}] cone {len(conewires)} wires, {len(patterns)} patterns ({time.time()-t0:.0f}s)", flush=True)
    allbits = sorted(set().union(*[set(p) for p in patterns])) if patterns else []
    hitsets = []
    CH = 4096
    for p in PRIMES:
        invmod = lambda z: pow(int(z) % p, p-2, p)
        bm = {v: int(bestval[v] % p) for v in range(NVARS)}
        rg = {v: (int(c % p), (invmod(c) if c % p else 0), [(int(cc % p), m) for cc, m in ts])
              for v, (c, ts) in exgate.items()}
        rl = {v: (bit, invmod(cbx), [(int(cc % p), m) for cc, m in lt])
              for v, (bit, cbx, lt) in exload.items()}
        hits = []
        for start in range(0, len(patterns), CH):
            batch = patterns[start:start+CH]
            L = len(batch)
            bitval = {}
            for b in allbits:
                col = np.zeros(L, dtype=np.int64)
                for i, pat in enumerate(batch):
                    if b in pat: col[i] = 1
                bitval[b] = col
            val = {}
            def getv(x):
                a = val.get(x)
                if a is not None: return a
                a = bitval.get(x)
                if a is not None: return a
                return bm[x]
            for v in conewires:
                if v in rl:
                    bit, invcbx, lt = rl[v]
                    rs = np.zeros(L, dtype=np.int64)
                    for c, m in lt:
                        t = np.full(L, c, dtype=np.int64)
                        for x in m: t = (t * (1 if x == bit else getv(x))) % p
                        rs = (rs + t) % p
                    lv = ((-rs) * invcbx) % p
                    ba = bitval.get(bit)
                    if ba is None:
                        bc = bm[bit]; val[v] = (lv * bc) % p if bc else np.zeros(L, dtype=np.int64)
                    else:
                        val[v] = (lv * ba) % p
                else:
                    coefp, invc, terms = rg[v]
                    if coefp == 0: continue
                    rs = np.zeros(L, dtype=np.int64)
                    for c, m in terms:
                        t = np.full(L, c, dtype=np.int64)
                        for x in m: t = (t * getv(x)) % p
                        rs = (rs + t) % p
                    val[v] = ((-rs) * invc) % p
            ok = np.ones(L, dtype=bool)
            for a in CHECKS:
                s = np.zeros(L, dtype=np.int64)
                for m, c in checkpoly[a]:
                    t = np.full(L, c % p, dtype=np.int64)
                    for x in m: t = (t * getv(x)) % p
                    s = (s + t) % p
                ok &= (s == 0)
            for idx in np.nonzero(ok)[0]:
                hits.append(batch[int(idx)])
        print(f"[{tag}] prime {p}: {len(hits)} hits ({time.time()-t0:.0f}s)", flush=True)
        hitsets.append([tuple(sorted(h)) for h in hits])
        if not hits: break
    if len(hitsets) < 2:
        common = []
    else:
        common = sorted(set(hitsets[0]) & set(hitsets[1]))
    print(f"[{tag}] common hits: {len(common)}: {common[:20]}", flush=True)
    json.dump([list(h) for h in common], open(f'hits_{tag}.json', 'w'))
    return common

def main():
    control = json.load(open('control_bits.json'))
    phase = sys.argv[1] if len(sys.argv) > 1 else 'allpairs'
    if phase == 'allpairs':
        pats = list({frozenset(c) for c in itertools.combinations(control, 2)})
        search(pats, "allpairs")
    elif phase == 'trip2in22':
        pats = set()
        for a, b in itertools.combinations(BITS22, 2):
            for c in control:
                if c != a and c != b: pats.add(frozenset([a, b, c]))
        search(list(pats), "trip2in22")
    elif phase == 'trip1in22':
        # 1 from 22, 2 from outside (activates a product of two pair-only bits, gated by a linear bit)
        out = [b for b in control if b not in set(BITS22)]
        pats = set()
        for a in BITS22:
            for b, c in itertools.combinations(out, 2):
                pats.add(frozenset([a, b, c]))
        print("trip1in22 count:", len(pats))
        search(list(pats), "trip1in22")

if __name__ == '__main__':
    main()
