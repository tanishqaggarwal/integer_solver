#!/usr/bin/env python3
"""Vectorized 2^22 enumeration over the 22 twist-bits (confluent evaluator).

Correct load-aware injection evaluator, reduced to the 22-bit forward cone
(~1865 wires). Evaluate all 2^22 patterns in numpy chunks mod a 31-bit prime;
keep patterns zeroing both check conditions (x_18274=x_9770, x_17728=x_3183).
Screen with two primes; verify survivors exactly in Z."""
import json, time, sys
import numpy as np
from collections import deque, defaultdict
from confluent_eval4 import build, boolean_vars
from propagate import load_atoms, atom_vars, NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
CHECKS = [1817, 30378, 44271]
PRIMES = [2147483647, 2147483629]
CHUNK = 1 << 15

def main():
    t0 = time.time()
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

    # cone
    deps = defaultdict(set)
    for v in gate:
        for c, m in gate[v][1]: deps[v].update(m)
        if v in loadinj:
            bit, cbx, lt = loadinj[v]; deps[v].add(bit)
            for c, m in lt: deps[v].update(m)
        deps[v].discard(v)
    users = defaultdict(list)
    for v in gate:
        for x in deps[v]: users[x].append(v)
    cone = set(BITS22); dq = deque(BITS22)
    while dq:
        x = dq.popleft()
        for u in users[x]:
            if u not in cone: cone.add(u); dq.append(u)
    conewires = [v for v in seq if v in cone and v in gate]
    print(f"cone {len(cone)}, wires {len(conewires)} ({time.time()-t0:.0f}s)", flush=True)

    # re-detect huge-atoms exactly (with atom index)
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

    # exact per-wire structures
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
    bitidx = {b: i for i, b in enumerate(BITS22)}
    N = 1 << 22
    hitsets = []

    for p in PRIMES:
        invmod = lambda z: pow(int(z) % p, p-2, p)
        bm = {v: int(bestval[v] % p) for v in range(NVARS)}
        rg = {v: (int(c % p), (invmod(c) if c % p else 0),
                  [(int(cc % p), m) for cc, m in ts]) for v, (c, ts) in exgate.items()}
        rl = {v: (bit, invmod(cbx), [(int(cc % p), m) for cc, m in lt])
              for v, (bit, cbx, lt) in exload.items()}
        hits = []
        for start in range(0, N, CHUNK):
            codes = np.arange(start, min(start+CHUNK, N), dtype=np.int64)
            L = len(codes)
            bitval = {b: ((codes >> bitidx[b]) & 1) for b in BITS22}
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
                        for x in m:
                            t = (t * (1 if x == bit else getv(x))) % p
                        rs = (rs + t) % p
                    lv = ((-rs) * invcbx) % p
                    ba = bitval.get(bit)
                    if ba is None:
                        bc = bm[bit]
                        val[v] = (lv * bc) % p if bc else np.zeros(L, dtype=np.int64)
                    else:
                        val[v] = (lv * ba) % p
                else:
                    coefp, invc, terms = rg[v]
                    if coefp == 0: continue
                    rs = np.zeros(L, dtype=np.int64)
                    for c, m in terms:
                        t = np.full(L, c, dtype=np.int64)
                        for x in m:
                            t = (t * getv(x)) % p
                        rs = (rs + t) % p
                    val[v] = ((-rs) * invc) % p
            ok = np.ones(L, dtype=bool)
            for a in CHECKS:
                s = np.zeros(L, dtype=np.int64)
                for m, c in checkpoly[a]:
                    t = np.full(L, c % p, dtype=np.int64)
                    for x in m:
                        t = (t * getv(x)) % p
                    s = (s + t) % p
                ok &= (s == 0)
            for idx in np.nonzero(ok)[0]:
                hits.append(int(codes[idx]))
            if (start // CHUNK) % 32 == 0:
                print(f"  p={p} {start}/{N} hits={len(hits)} ({time.time()-t0:.0f}s)", flush=True)
        print(f"prime {p}: {len(hits)} hits ({time.time()-t0:.0f}s)", flush=True)
        hitsets.append(set(hits))
        if len(hits) == 0: break
    common = set.intersection(*hitsets) if hitsets and all(hitsets) else (hitsets[0] if hitsets else set())
    print(f"common hits: {len(common)}", flush=True)
    json.dump([[b for i, b in enumerate(BITS22) if (h >> i) & 1] for h in sorted(common)],
              open('confluent_enum_hits.json', 'w'))
    print(f"wrote confluent_enum_hits.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
