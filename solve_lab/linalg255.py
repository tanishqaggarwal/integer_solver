#!/usr/bin/env python3
"""CORRECTED linear-algebra attack over ALL 255 control bits (mod P, artifact-free).

Session-5's linalg fixed the 22 bits to 0 and found rank 233 => B=0, which is
CONTRADICTORY (the witness has the 22 bits != 0). Fix: treat all 255 control bits
as variables. Compute each atom's mod-P single-flip response, keep atoms that are
LINEAR in the bits (multi-config test), Gaussian-eliminate over GF(P). The witness
(255-bit 0/1 vector) must satisfy every linear check, so it lies in the affine
solution space V. Report dim(V) = 255 - rank. If small, enumerate 0/1 points in V,
apply the twist, verify in Z.
"""
import json, time, sys
import multiprocessing as mp
from confluent_eval5 import build5, make_forward
from propagate import NVARS
from modp import P, inv

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
TWIST = {1817, 30378, 40782, 44271}

_G = {}
def init():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    _G.update(A=A, solveP=solveP, bm=bm, kind=kind, info=info, seq=seq, bestval=bestval)

def atomvals(setbits):
    val = _G['solveP'](list(_G['bm']), list(setbits))
    A = _G['A']; out = []
    for poly in A:
        s = 0
        for m, c in poly.items():
            t = c % P
            for x in m: t = (t*val[x]) % P
            s = (s+t) % P
        out.append(s)
    return out

def worker(bit):
    return bit, atomvals([] if bit is None else [bit])

def main():
    t0 = time.time()
    init()
    control = json.load(open('control_bits.json'))
    print(f"control bits: {len(control)}", flush=True)
    res = {}
    with mp.Pool(6, initializer=init) as pool:
        for k, (bit, av) in enumerate(pool.imap_unordered(worker, [None] + control)):
            res[bit] = av
            if (k+1) % 40 == 0: print(f"  {k+1}/{len(control)+1} responses ({time.time()-t0:.0f}s)", flush=True)
    base = res[None]; NA = len(base)
    col = {b: i for i, b in enumerate(control)}; n = len(control)
    # linearity filter: multi-bit config, mod P
    testcfg = control[::2]
    av_test = atomvals(testcfg)
    islin = []
    for a in range(NA):
        pred = base[a]
        for b in testcfg: pred = (pred + res[b][a] - base[a]) % P
        islin.append(pred % P == av_test[a] % P)
    # second config to reduce false "linear"
    testcfg2 = control[1::3]
    av_test2 = atomvals(testcfg2)
    for a in range(NA):
        if not islin[a]: continue
        pred = base[a]
        for b in testcfg2: pred = (pred + res[b][a] - base[a]) % P
        if pred % P != av_test2[a] % P: islin[a] = False
    print(f"atoms linear in the 255 bits: {sum(islin)}/{NA} ({time.time()-t0:.0f}s)", flush=True)

    # Build linear system A_lin . x = rhs over GF(P): sum_b delta_{a,b} x_b = -base[a]
    # Gaussian elimination tracking rank; also detect inconsistency at x that must be 0/1.
    pivots = {}   # col -> (row dict, rhs)
    rowsused = 0; incons = 0
    twist_rows = 0
    for a in range(NA):
        if not islin[a]: continue
        row = {}
        for b in control:
            d = (res[b][a] - base[a]) % P
            if d: row[col[b]] = d
        rhs = (-base[a]) % P
        if not row and rhs == 0: continue
        rowsused += 1
        if a in TWIST: twist_rows += 1
        r = dict(row); rr = rhs
        while True:
            pc = [c for c in r if c in pivots]
            if not pc: break
            c = pc[0]; f = r[c]; prow, prhs = pivots[c]
            for k, v in prow.items():
                nv = (r.get(k, 0) - f*v) % P
                if nv: r[k] = nv
                elif k in r: del r[k]
            rr = (rr - f*prhs) % P
        if not r:
            if rr != 0: incons += 1
            continue
        pcn = min(r); fi = inv(r[pcn])
        pivots[pcn] = ({k: (v*fi) % P for k, v in r.items()}, (rr*fi) % P)
    rank = len(pivots)
    print(f"linear rows used {rowsused} (twist among them: {twist_rows}); RANK {rank}/{n}; free dim {n-rank}; inconsistent {incons} ({time.time()-t0:.0f}s)", flush=True)
    free = [b for i, b in enumerate(control) if i not in pivots]
    print(f"free bit-columns ({len(free)}): {sorted(free)}", flush=True)
    # Is BITS22 exactly the free set?
    print(f"free set == BITS22 ? {set(free) == set(BITS22)}; overlap {len(set(free)&set(BITS22))}/{len(BITS22)}", flush=True)
    # particular solution (free = 0) then check 0/1
    sol = [0]*n
    for pc, (prow, prhs) in pivots.items():
        val = prhs
        for k, v in prow.items():
            if k != pc: val = (val - v*sol[k]) % P
        sol[pc] = val % P
    nonbin = sum(1 for x in sol if x not in (0,1))
    print(f"particular solution (free=0): {sum(1 for x in sol if x==1)} ones, {nonbin} non-binary", flush=True)
    json.dump({'rank': rank, 'free': sorted(free), 'nfree': n-rank,
               'freeIsBITS22': set(free)==set(BITS22)}, open('linalg255.json','w'))
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
