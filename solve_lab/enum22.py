#!/usr/bin/env python3
"""Solve by reducing to the 22 nonlinear bits.

The 233 x_18274-side bits are slaved to the 22 bits: B(A) = -M^{-1} base(A), where
M is the (full-rank) 232-bit response matrix of 233 independent linear atoms and
base(A) is their value at (22=A, 232=0). B(A) is 0/1 only for the true A. Enumerate
22-bit patterns A (small subsets first), compute B(A), filter by 0/1, verify in Z."""
import json, time, sys, itertools
import multiprocessing as mp
from confluent_eval5 import build5, make_forward
from propagate import NVARS
from modp import P, inv

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

_G = {}
def init():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    bm = [x % P for x in bestval]
    _G.update(A=A, kind=kind, info=info, seq=seq, bm=bm, bestval=bestval)

def fwd(setbits):
    kind = _G['kind']; info = _G['info']; seq = _G['seq']
    val = list(_G['bm'])
    for b in setbits: val[b] = 1
    for v in seq:
        k = kind[v]
        if k == 'gate':
            coef, terms = info[v]
            if coef % P == 0: continue
            rs = 0
            for c, m in terms:
                t = c % P
                for x in m: t = (t*val[x]) % P
                rs = (rs+t) % P
            val[v] = (-rs * inv(coef)) % P
        elif k == 'load':
            bit, cbx, lt = info[v]
            if val[bit] == 0: val[v] = 0
            else:
                rs = 0
                for c, m in lt:
                    t = c % P
                    for x in m: t = (t*(1 if x == bit else val[x])) % P
                    rs = (rs+t) % P
                val[v] = (-rs * inv((cbx*val[bit]) % P)) % P
        elif k == 'div':
            c, u, rest = info[v]
            if val[u] == 0: val[v] = 0; continue
            rs = 0
            for cc, m in rest:
                t = cc % P
                for x in m: t = (t*val[x]) % P
                rs = (rs+t) % P
            val[v] = (-rs * inv((c*val[u]) % P)) % P
    return val

def atomvals(val, idxs=None):
    A = _G['A']
    if idxs is None: idxs = range(len(A))
    out = []
    for a in idxs:
        s = 0
        for m, c in A[a].items():
            t = c % P
            for x in m: t = (t*val[x]) % P
            s = (s+t) % P
        out.append(s)
    return out

def resp_worker(bit):
    return bit, atomvals(fwd([bit] if bit is not None else []))

def base_worker(args):
    Aset, pivot_atoms = args
    val = fwd(list(Aset))
    return tuple(Aset), atomvals(val, pivot_atoms)

def main():
    t0 = time.time()
    init()
    control = json.load(open('control_bits.json'))
    bits232 = [b for b in control if b not in set(BITS22)]
    # 1) responses of all atoms to the 232 bits
    res = {}
    with mp.Pool(6, initializer=init) as pool:
        for bit, av in pool.imap_unordered(resp_worker, [None] + bits232):
            res[bit] = av
    base0 = res[None]; NA = len(base0)
    # linearity filter
    testcfg = bits232[::2]; av_test = atomvals(fwd(testcfg))
    islin = []
    for a in range(NA):
        pred = base0[a]
        for b in testcfg: pred = (pred + res[b][a] - base0[a]) % P
        islin.append(pred % P == av_test[a] % P)
    print(f"linear atoms {sum(islin)}/{NA} ({time.time()-t0:.0f}s)", flush=True)
    # 2) build M from 233 independent linear atoms (Gaussian, track pivot atoms)
    col = {b: i for i, b in enumerate(bits232)}; n = len(bits232)
    pivrow = {}      # pivot col -> (reduced row dict, atom index)
    pivot_atoms = []
    for a in range(NA):
        if not islin[a]: continue
        row = {}
        for b in bits232:
            d = (res[b][a] - base0[a]) % P
            if d: row[col[b]] = d
        if not row: continue
        r = dict(row)
        while True:
            pc = [c for c in r if c in pivrow]
            if not pc: break
            c = pc[0]; f = r[c]; prow, _ = pivrow[c]
            for k, v in prow.items():
                nv = (r.get(k, 0) - f*v) % P
                if nv: r[k] = nv
                elif k in r: del r[k]
        if not r: continue
        pcnew = min(r); fi = inv(r[pcnew])
        pivrow[pcnew] = ({k: (v*fi) % P for k, v in r.items()}, a)
        pivot_atoms.append(a)
        if len(pivot_atoms) == n: break
    print(f"pivot atoms {len(pivot_atoms)} (need {n}) ({time.time()-t0:.0f}s)", flush=True)
    if len(pivot_atoms) < n:
        print("system under-determined; cannot slave 232-part uniquely", flush=True); return
    # M = pivot atoms' 232-response rows (n x n), invert mod P
    import numpy as np
    M = [[ (res[bits232[j]][pivot_atoms[i]] - base0[pivot_atoms[i]]) % P for j in range(n)] for i in range(n)]
    # invert M mod P via Gauss-Jordan
    Minv = gj_inverse(M, P)
    if Minv is None:
        print("M singular", flush=True); return
    print(f"M inverted ({time.time()-t0:.0f}s)", flush=True)

    def solveB(baseA):   # baseA: list of pivot-atom values at (22=A,232=0); B = -Minv @ baseA
        B = [0]*n
        for i in range(n):
            s = 0; row = Minv[i]
            for j in range(n): s = (s + row[j]*baseA[j]) % P
            B[i] = (-s) % P
        return B

    # 3) enumerate small 22-subsets
    subsets = [()]
    for k in range(1, 4):
        subsets += list(itertools.combinations(BITS22, k))
    print(f"testing {len(subsets)} small 22-subsets ({time.time()-t0:.0f}s)", flush=True)
    tasks = [(s, pivot_atoms) for s in subsets]
    hits = []
    with mp.Pool(6, initializer=init) as pool:
        for k, (Aset, baseA) in enumerate(pool.imap_unordered(base_worker, tasks)):
            B = solveB(baseA)
            nonbin = sum(1 for x in B if x not in (0, 1))
            if nonbin == 0:
                ones232 = [bits232[i] for i in range(n) if B[i] == 1]
                allones = list(Aset) + ones232
                hits.append(allones)
                print(f"  0/1 B for A={Aset}: {len(ones232)} 232-bits set", flush=True)
            if (k+1) % 200 == 0:
                print(f"  ...{k+1}/{len(subsets)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"0/1 candidates: {len(hits)} ({time.time()-t0:.0f}s)", flush=True)
    # verify in Z
    solve = make_forward(_G['kind'], _G['info'], _G['seq'], _G['bestval'])
    for allones in hits:
        vz = solve(list(_G['bestval']), allones)
        vio = 0
        for poly in _G['A']:
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= vz[x]
                s += t
            if s: vio += 1
        print(f"  Z-verify {sorted(allones)[:8]}...: {vio} violated", flush=True)
        if vio == 0:
            json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json', 'w'))
            print("  *** SOLVED! wrote cand_SOLVED.json ***", flush=True); return

def gj_inverse(M, P):
    n = len(M)
    A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] % P), None)
        if piv is None: return None
        A[c], A[piv] = A[piv], A[c]
        f = inv(A[c][c])
        A[c] = [(x*f) % P for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                g = A[r][c]
                A[r] = [(A[r][k] - g*A[c][k]) % P for k in range(2*n)]
    return [row[n:] for row in A]

if __name__ == '__main__':
    main()
