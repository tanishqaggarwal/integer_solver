#!/usr/bin/env python3
"""Tabulate the 22-side: for ALL 2^22 patterns over BITS22 compute
(x_9770, x_3183) mod two 31-bit primes, using a cone-restricted vectorized eval.
Check collision with the B=0 target (x_18274(0), x_17728(0)) mod the primes.
Save the table (packed) for later meet-in-the-middle vs the B-side.
"""
import json, time, sys
import numpy as np
from collections import deque
from confluent_eval5 import build5, make_forward
from propagate import atom_vars, NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
PRIMES = [2147483647, 2147483629]
TARGETVARS = [9770, 3183]

def prep():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    deps = {}
    for v in defset:
        d = set(); k = kind[v]
        if k == 'gate':
            for c, m in info[v][1]: d.update(m)
        elif k == 'load':
            bit, cbx, lt = info[v]; d.add(bit)
            for c, m in lt: d.update(m)
        elif k == 'div':
            c, u, rest = info[v]; d.add(u)
            for cc, m in rest: d.update(m)
        d.discard(v); deps[v] = d
    seen = set(TARGETVARS); dq = deque(TARGETVARS)
    while dq:
        v = dq.popleft()
        for x in deps.get(v, ()):
            if x not in seen: seen.add(x); dq.append(x)
    conewires = [v for v in seq if v in seen and v in defset]
    # B=0 target from full Z eval
    solveZ = make_forward(kind, info, seq, bestval)
    base = solveZ(list(bestval), [])
    tgt = {18274: base[18274], 17728: base[17728]}
    return A, kind, info, conewires, bestval, seen, tgt

def eval_batch(batch, conewires, kind, info, bm, p):
    L = len(batch)
    bitcols = {}
    for i, mask in enumerate(batch):
        pass
    # build bit indicator columns
    bitval = {}
    for bi, b in enumerate(BITS22):
        col = ((batch >> bi) & 1).astype(np.int64)
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
            ba = bitval.get(bit); invc = pow(int(cbx), p-2, p); lv = ((-rs) * invc) % p
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
            if np.isscalar(uv):
                den = (int(c) % p) * (uv % p) % p
                val[v] = ((-rs) * (pow(int(den), p-2, p) if den else 0)) % p
            else:
                den = (int(c) % p * uv) % p
                iv = vinv(den, p); iv = np.where(den == 0, 0, iv)
                val[v] = ((-rs) * iv) % p
    return getv(9770) % p, getv(3183) % p

def vinv(arr, p):
    r = np.ones_like(arr); b = arr % p; e = p - 2
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r

def main():
    t0 = time.time()
    A, kind, info, conewires, bestval, seen, tgt = prep()
    ndiv = sum(1 for v in conewires if kind[v] == 'div')
    print(f"22-side cone: {len(seen)} vars, {len(conewires)} wires ({ndiv} div) ({time.time()-t0:.0f}s)", flush=True)
    N = 1 << 22
    CH = 1 << 16
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else N
    for p in PRIMES:
        bm = {v: int(bestval[v] % p) for v in range(NVARS)}
        t1 = tgt[18274] % p; t2 = tgt[17728] % p
        arr9770 = np.empty(min(limit, N), dtype=np.int64)
        arr3183 = np.empty(min(limit, N), dtype=np.int64)
        nhit = 0
        for start in range(0, min(limit, N), CH):
            end = min(start + CH, min(limit, N))
            batch = np.arange(start, end, dtype=np.int64)
            v9770, v3183 = eval_batch(batch, conewires, kind, info, bm, p)
            arr9770[start:end] = v9770; arr3183[start:end] = v3183
            hit = np.nonzero((v9770 == t1) & (v3183 == t2))[0]
            nhit += len(hit)
            if len(hit):
                for h in hit: print(f"  [{p}] B=0 collision at pattern {start+int(h):#x}", flush=True)
            if start % (CH*16) == 0:
                print(f"  [{p}] {end}/{min(limit,N)} ({time.time()-t0:.0f}s)", flush=True)
        # also report how many match on x_9770 alone (coordinate 1)
        n1 = int(np.count_nonzero(arr9770 == t1)); n2 = int(np.count_nonzero(arr3183 == t2))
        print(f"[{p}] target t1={t1} t2={t2}: single-coord matches x9770={n1} x3183={n2}, BOTH={nhit} ({time.time()-t0:.0f}s)", flush=True)
        np.save(f'tab22_9770_{p}.npy', arr9770); np.save(f'tab22_3183_{p}.npy', arr3183)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
