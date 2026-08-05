#!/usr/bin/env python3
"""Vectorized bit-pattern search with the CORRECTED v5 evaluator (division wires).

v4 froze x_18274; v5 computes it (and 1374 other division wires), so x_18274 is
bit-movable and the earlier 0-hits are void. This re-runs the search with the
correct model. Handles gate / load / div wires vectorized (div uses elementwise
modpow inverse). Two 31-bit primes; verify hits exactly in Z."""
import json, time, sys, itertools
import numpy as np
from collections import defaultdict, deque
from confluent_eval5 import build5
from propagate import load_atoms, atom_vars, NVARS

CHECKS = [1817, 30378, 44271]
PRIMES = [2147483647, 2147483629]
BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def vinv(arr, p):
    r = np.ones_like(arr); b = arr % p; e = p - 2
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r

def prep():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    # deps for cone
    deps = {}
    for v in defset:
        d = set()
        k = kind[v]
        if k == 'gate':
            for c, m in info[v][1]: d.update(m)
        elif k == 'load':
            bit, cbx, lt = info[v]; d.add(bit)
            for c, m in lt: d.update(m)
        elif k == 'div':
            c, u, rest = info[v]; d.add(u)
            for cc, m in rest: d.update(m)
        d.discard(v); deps[v] = d
    checkvars = set()
    for a in CHECKS: checkvars.update(atom_vars(A[a]))
    seen = set(checkvars); dq = deque(checkvars)
    while dq:
        v = dq.popleft()
        for x in deps.get(v, ()):
            if x not in seen: seen.add(x); dq.append(x)
    conewires = [v for v in seq if v in seen and v in defset]
    checkpoly = {a: list(A[a].items()) for a in CHECKS}
    return A, kind, info, conewires, checkpoly, bestval, len(seen)

def search(patterns, tag):
    t0 = time.time()
    A, kind, info, conewires, checkpoly, bestval, conesz = prep()
    ndiv = sum(1 for v in conewires if kind[v] == 'div')
    print(f"[{tag}] cone {conesz}, wires {len(conewires)} ({ndiv} div), {len(patterns)} patterns ({time.time()-t0:.0f}s)", flush=True)
    allbits = sorted(set().union(*[set(p) for p in patterns])) if patterns else []
    hitsets = []
    CH = 2048
    for p in PRIMES:
        bm = {v: int(bestval[v] % p) for v in range(NVARS)}
        hits = []
        for start in range(0, len(patterns), CH):
            batch = patterns[start:start+CH]; L = len(batch)
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
                k = kind[v]
                if k == 'gate':
                    coef, terms = info[v]
                    if coef % p == 0: continue
                    rs = np.zeros(L, dtype=np.int64)
                    for c, m in terms:
                        t = np.full(L, c % p, dtype=np.int64)
                        for x in m: t = (t * getv(x)) % p
                        rs = (rs + t) % p
                    val[v] = ((-rs) * pow(int(coef), p-2, p)) % p
                elif k == 'load':
                    bit, cbx, lt = info[v]
                    rs = np.zeros(L, dtype=np.int64)
                    for c, m in lt:
                        t = np.full(L, c % p, dtype=np.int64)
                        for x in m: t = (t * (1 if x == bit else getv(x))) % p
                        rs = (rs + t) % p
                    ba = bitval.get(bit)
                    invc = pow(int(cbx), p-2, p)
                    lv = ((-rs) * invc) % p
                    if ba is None:
                        bc = bm[bit]; val[v] = (lv * bc) % p if bc else np.zeros(L, dtype=np.int64)
                    else:
                        val[v] = (lv * ba) % p
                elif k == 'div':
                    c, u, rest = info[v]
                    rs = np.zeros(L, dtype=np.int64)
                    for cc, m in rest:
                        t = np.full(L, cc % p, dtype=np.int64)
                        for x in m: t = (t * getv(x)) % p
                        rs = (rs + t) % p
                    uv = getv(u)
                    den = (int(c) % p) * (uv % p) % p if np.isscalar(uv) else (int(c) % p * uv) % p
                    if np.isscalar(den):
                        val[v] = ((-rs) * (pow(int(den), p-2, p) if den else 0)) % p
                    else:
                        iv = vinv(den, p)
                        iv = np.where(den == 0, 0, iv)
                        val[v] = ((-rs) * iv) % p
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
        hitsets.append(set(tuple(sorted(h)) for h in hits))
        if not hits: break
    common = sorted(set.intersection(*hitsets)) if len(hitsets) == 2 else []
    print(f"[{tag}] common hits: {len(common)}: {common[:20]}", flush=True)
    json.dump([list(h) for h in common], open(f'hits5_{tag}.json', 'w'))
    return common

def main():
    control = json.load(open('control_bits.json'))
    phase = sys.argv[1] if len(sys.argv) > 1 else 'allpairs'
    if phase == 'singles':
        search([frozenset([b]) for b in control], 'singles')
    elif phase == 'allpairs':
        search(list({frozenset(c) for c in itertools.combinations(control, 2)}), 'allpairs')

if __name__ == '__main__':
    main()
