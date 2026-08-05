#!/usr/bin/env python3
"""Test the IMAGE SIZE of the 233-side: sample many random B and count distinct
(x_18274, x_17728) mod p. The 22-side collapsed to 285 pairs; if the 233-side also
has a small image, the whole match is a tiny finite intersection."""
import json, time
import numpy as np
from mitm_lowB import prep, vinv
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
p = 2147483647

def main():
    t0 = time.time()
    A, kind, info, conewires, bestval = prep([18274, 17728])
    bm = {v: int(bestval[v] % p) for v in range(NVARS)}
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]
    print(f"233-cone {len(conewires)} wires ({time.time()-t0:.0f}s)", flush=True)

    def eval_batch(bmasks):
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

    # sample random B of various weights (deterministic LCG)
    st = 999
    def rnd():
        nonlocal st; st = (st*1103515245+12345) & 0x7fffffff; return st
    seen18 = set(); seen17 = set(); seenpair = set()
    NB = 40000
    B = []
    for _ in range(NB):
        w = 1 + rnd() % 60
        s = set()
        for _ in range(w): s.add(bits233[rnd() % len(bits233)])
        B.append(frozenset(s))
    CH = 20000
    for start in range(0, len(B), CH):
        batch = B[start:start+CH]
        x18, x17 = eval_batch(batch)
        for i in range(len(batch)):
            seen18.add(int(x18[i])); seen17.add(int(x17[i]))
            seenpair.add((int(x18[i]), int(x17[i])))
    print(f"sampled {NB} random B (weight 1-60):", flush=True)
    print(f"  distinct x_18274: {len(seen18)}", flush=True)
    print(f"  distinct x_17728: {len(seen17)}", flush=True)
    print(f"  distinct (x_18274,x_17728) pairs: {len(seenpair)}", flush=True)
    # compare with 22-side achievable set
    a9 = np.load('tab22_9770_2147483647.npy'); a3 = np.load('tab22_3183_2147483647.npy')
    pairs22 = set(zip(a9.tolist(), a3.tolist()))
    inter = seenpair & pairs22
    print(f"  22-side pairs: {len(pairs22)}; INTERSECTION with sampled 233-pairs: {len(inter)}", flush=True)
    if inter:
        print(f"  *** intersecting pairs: {list(inter)[:5]} ***", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
