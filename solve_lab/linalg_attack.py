#!/usr/bin/env python3
"""Linear-algebra attack over GF(P).

For a bit-config restricted to the 232 x_18274-side bits, every atom that is
linear in those bits satisfies  base_a + sum_i delta_{a,i} b_i = 0.  base_a = 0
for all atoms except the twist (at all-0). Collect all such linear constraints
and Gaussian-eliminate over GF(P): if full rank the unique solution is the
witness's 232-part (0/1 read mod P). Verify exactly with v5."""
import json, time, sys
import multiprocessing as mp
from confluent_eval5 import build5
from propagate import NVARS
from modp import P, inv

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
    A = _G['A']; kind = _G['kind']; info = _G['info']; seq = _G['seq']
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

def atomvals(val):
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
    return bit, atomvals(fwd([bit] if bit is not None else []))

def main():
    t0 = time.time()
    init()
    control = json.load(open('control_bits.json'))
    BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
    bits232 = [b for b in control if b not in BITS22]
    print(f"solving over {len(bits232)} bits", flush=True)
    tasks = [None] + bits232
    res = {}
    with mp.Pool(6, initializer=init) as pool:
        for k, (bit, av) in enumerate(pool.imap_unordered(worker, tasks)):
            res[bit] = av
            if (k+1) % 40 == 0: print(f"  {k+1}/{len(tasks)} responses ({time.time()-t0:.0f}s)", flush=True)
    base = res[None]
    NA = len(base)
    # linearity test: compare actual value at a multi-bit config to the linear prediction
    testcfg = bits232[::2]   # every other bit
    av_test = atomvals(fwd(testcfg))
    testset = set(testcfg)
    is_linear = [False]*NA
    for a in range(NA):
        pred = base[a]
        for b in testcfg: pred = (pred + res[b][a] - base[a]) % P
        is_linear[a] = (pred % P == av_test[a] % P)
    nlin = sum(is_linear)
    print(f"atoms linear in the 232 bits: {nlin}/{NA} ({time.time()-t0:.0f}s)", flush=True)
    # build sparse linear system from LINEAR atoms only
    col = {b: i for i, b in enumerate(bits232)}
    rows = []   # (dict col->coef, rhs)
    for a in range(NA):
        if not is_linear[a]: continue
        row = {}
        for b in bits232:
            d = (res[b][a] - base[a]) % P
            if d: row[col[b]] = d
        if row or base[a]:
            rows.append((row, (-base[a]) % P))
    print(f"linear constraint rows: {len(rows)} ({time.time()-t0:.0f}s)", flush=True)

    # Gaussian elimination over GF(P)
    n = len(bits232)
    pivots = {}   # col -> row (as dict), rhs
    used = []
    incons = 0
    for row, rhs in rows:
        row = dict(row); r = rhs
        # reduce by existing pivots (loop until no pivot col remains)
        while True:
            pc_in = [c for c in row if c in pivots]
            if not pc_in: break
            c = pc_in[0]; f = row[c]
            prow, prhs = pivots[c]
            for k, v in prow.items():
                nv = (row.get(k, 0) - f*v) % P
                if nv: row[k] = nv
                elif k in row: del row[k]
            r = (r - f*prhs) % P
        if not row:
            if r != 0: incons += 1
            continue
        # pick a pivot col
        pc = min(row)
        f = inv(row[pc])
        prow = {k: (v*f) % P for k, v in row.items()}
        pivots[pc] = (prow, (r*f) % P)
    rank = len(pivots)
    print(f"rank {rank} / {n} vars, inconsistent rows {incons} ({time.time()-t0:.0f}s)", flush=True)
    # back-substitute: free vars = 0
    sol = [0]*n
    for pc, (prow, prhs) in pivots.items():
        val = prhs
        for k, v in prow.items():
            if k != pc: val = (val - v*sol[k]) % P
        sol[pc] = val % P
    ones = [bits232[i] for i in range(n) if sol[i] == 1]
    nonbin = [bits232[i] for i in range(n) if sol[i] not in (0, 1)]
    print(f"solution: {len(ones)} bits=1, {len(nonbin)} non-binary", flush=True)
    print(f"  bits=1: {sorted(ones)[:40]}", flush=True)
    json.dump({'ones': sorted(ones), 'nonbin_count': len(nonbin), 'rank': rank}, open('linalg_sol.json', 'w'))
    # verify if binary
    if not nonbin:
        val = fwd(ones)
        av = atomvals(val)
        viol = sum(1 for x in av if x != 0)
        print(f"  mod-P violated atoms: {viol}", flush=True)
        if viol == 0:
            from confluent_eval5 import make_forward
            solve = make_forward(_G['kind'], _G['info'], _G['seq'], _G['bestval'])
            vz = solve(list(_G['bestval']), ones)
            vio = 0
            for poly in _G['A']:
                s = 0
                for m, c in poly.items():
                    t = c
                    for x in m: t *= vz[x]
                    s += t
                if s: vio += 1
            print(f"  *** Z violated atoms: {vio} ***", flush=True)
            if vio == 0:
                json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json', 'w'))
                print("  *** SOLVED! ***", flush=True)

if __name__ == '__main__':
    main()
