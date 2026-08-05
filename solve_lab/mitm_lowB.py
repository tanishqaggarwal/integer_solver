#!/usr/bin/env python3
"""MITM: hash the full 2^22 22-side table {(x_9770(A),x_3183(A)) mod p -> A}, then
sweep low-weight 233-side B, computing (x_18274(B),x_17728(B)) mod p via a vectorized
233-cone eval, and look up. A hit => x_18274(B)=x_9770(A) & x_17728(B)=x_3183(A) for
SOME A (any weight) with a sparse B. Z-verify hits. More powerful than prior symmetric
(total-weight<=3) searches for asymmetric witnesses."""
import json, time, sys, itertools
import numpy as np
from collections import deque
from confluent_eval5 import build5, make_forward
from propagate import atom_vars, NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
p = 2147483647   # matches tab22_*_2147483647.npy

def prep(targets):
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
    seen = set(targets); dq = deque(targets)
    while dq:
        v = dq.popleft()
        for x in deps.get(v, ()):
            if x not in seen: seen.add(x); dq.append(x)
    conewires = [v for v in seq if v in seen and v in defset]
    return A, kind, info, conewires, bestval

def vinv(arr, pp):
    r = np.ones_like(arr); b = arr % pp; e = pp - 2
    while e:
        if e & 1: r = (r*b) % pp
        b = (b*b) % pp; e >>= 1
    return r

def main():
    t0 = time.time()
    a9770 = np.load('tab22_9770_2147483647.npy'); a3183 = np.load('tab22_3183_2147483647.npy')
    print(f"loaded 22-side table {len(a9770)} ({time.time()-t0:.0f}s)", flush=True)
    # hash: key = x9770*p + x3183  (fits in int64? p~2^31, product ~2^62 ok)
    keys = (a9770.astype(np.int64) * p + a3183.astype(np.int64))
    hashset = {}
    for A in range(len(keys)):
        hashset[int(keys[A])] = A
    print(f"built hash {len(hashset)} ({time.time()-t0:.0f}s)", flush=True)

    A, kind, info, conewires, bestval = prep([18274, 17728])
    bm = {v: int(bestval[v] % p) for v in range(NVARS)}
    ndiv = sum(1 for v in conewires if kind[v] == 'div')
    print(f"233-cone {len(conewires)} wires ({ndiv} div) ({time.time()-t0:.0f}s)", flush=True)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]

    def eval_batch(bmasks, blist):
        # bmasks: list of frozensets of set bits; compute x18274,x17728 mod p vectorized
        L = len(bmasks)
        allb = sorted(set().union(*bmasks)) if bmasks else []
        bitval = {}
        for b in allb:
            col = np.zeros(L, dtype=np.int64)
            for i, ms in enumerate(bmasks):
                if b in ms: col[i] = 1
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
                    for x in m: t = (t*getv(x)) % p
                    rs = (rs+t) % p
                val[v] = ((-rs) * pow(int(coef), p-2, p)) % p
            elif k == 'load':
                bit, cbx, lt = info[v]
                rs = np.zeros(L, dtype=np.int64)
                for c, m in lt:
                    t = np.full(L, c % p, dtype=np.int64)
                    for x in m: t = (t*(1 if x == bit else getv(x))) % p
                    rs = (rs+t) % p
                ba = bitval.get(bit); invc = pow(int(cbx), p-2, p); lv = ((-rs)*invc) % p
                if ba is None:
                    bc = bm[bit]; val[v] = (lv*bc) % p if bc else np.zeros(L, dtype=np.int64)
                else: val[v] = (lv*ba) % p
            elif k == 'div':
                c, u, rest = info[v]
                rs = np.zeros(L, dtype=np.int64)
                for cc, m in rest:
                    t = np.full(L, cc % p, dtype=np.int64)
                    for x in m: t = (t*getv(x)) % p
                    rs = (rs+t) % p
                uv = getv(u)
                if np.isscalar(uv):
                    den = (int(c) % p)*(uv % p) % p
                    val[v] = ((-rs)*(pow(int(den), p-2, p) if den else 0)) % p
                else:
                    den = (int(c) % p*uv) % p; iv = vinv(den, p); iv = np.where(den == 0, 0, iv)
                    val[v] = ((-rs)*iv) % p
        return getv(18274) % p, getv(17728) % p

    maxw = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    # enumerate low-weight B
    patterns = [frozenset()]
    for w in range(1, maxw+1):
        patterns += [frozenset(c) for c in itertools.combinations(bits233, w)]
    print(f"sweeping {len(patterns)} B-patterns (weight<= {maxw}) ({time.time()-t0:.0f}s)", flush=True)
    CH = 20000
    hits = []
    for start in range(0, len(patterns), CH):
        batch = patterns[start:start+CH]
        x18, x17 = eval_batch(batch, batch)
        k = (x18.astype(np.int64)*p + x17.astype(np.int64))
        for i in range(len(batch)):
            if x18[i] == 0 and x17[i] == 0: continue   # skip degenerate zero-collision
            A = hashset.get(int(k[i]))
            if A is not None:
                hits.append((batch[i], A))
                print(f"  NONZERO HIT B={sorted(batch[i])} <-> A-index {A} vals({x18[i]},{x17[i]}) ({time.time()-t0:.0f}s)", flush=True)
        if start % (CH*10) == 0:
            print(f"  ...{start+len(batch)}/{len(patterns)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"total hits: {len(hits)} ({time.time()-t0:.0f}s)", flush=True)
    # Z-verify hits
    if hits:
        solveZ = make_forward(kind, info, ([v for v in json.load(open('eval_order.json'))['order']]), bestval)
        # rebuild proper seq
        A2, kind2, info2, seq0, bv2, nc = build5()
        order = json.load(open('eval_order.json'))['order']
        defset = set(v for v in kind2 if kind2[v] != 'const')
        seq = [v for v in order if v in defset and v not in (9770, 3183)]
        seq += [v for v in (9770, 3183) if v in defset]
        seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
        solveZ = make_forward(kind2, info2, seq, bv2)
        A22bits = list(np.load('tab22_9770_2147483647.npy'))  # placeholder
        for Bset, Aidx in hits[:20]:
            # reconstruct A bits from index
            Abits = [BITS22[i] for i in range(22) if (Aidx >> i) & 1]
            allbits = list(Bset) + Abits
            vz = solveZ(list(bv2), allbits)
            vio = 0
            for poly in A2:
                s = 0
                for m, c in poly.items():
                    t = c
                    for x in m: t *= vz[x]
                    s += t
                if s: vio += 1
            print(f"  Z-verify A={Abits} B={sorted(Bset)}: {vio} violated", flush=True)
            if vio == 0:
                json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json','w'))
                print("  *** SOLVED! ***", flush=True); return
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
