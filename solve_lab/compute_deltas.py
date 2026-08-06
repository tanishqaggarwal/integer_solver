#!/usr/bin/env python3
"""Compute exact Z per-bit deltas of the two twist residuals (atoms 1817, 44271)
and check linearity. If (near-)linear, the twist is a 2-target subset-sum:
find bit set S with sum(delta1)= -R1(0) and sum(delta2) = -R2(0)."""
import json, time, sys
import multiprocessing as mp
from confluent_eval5 import build5, make_forward
from propagate import NVARS

_G = {}
def init():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    _G.update(A=A, solve=solve, bestval=bestval)

def resid(val, a):
    s = 0
    for m, c in _G['A'][a].items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

def evalset(setbits):
    val = _G['solve'](list(_G['bestval']), list(setbits))
    return resid(val, 1817), resid(val, 44271)

def worker(setbits):
    return setbits, evalset(setbits)

def main():
    init()
    control = json.load(open('control_bits.json'))
    base = evalset([])
    print(f"R1(0)={base[0]}", flush=True)
    print(f"R2(0)={base[1]}", flush=True)
    tasks = [tuple()] + [(b,) for b in control]
    d1 = {}; d2 = {}
    with mp.Pool(6, initializer=init) as pool:
        for sb, (r1, r2) in pool.imap_unordered(worker, tasks):
            if sb == tuple(): continue
            d1[sb[0]] = r1 - base[0]; d2[sb[0]] = r2 - base[1]
    movers1 = [b for b in control if d1[b] != 0]
    movers2 = [b for b in control if d2[b] != 0]
    print(f"bits moving R1: {len(movers1)}, R2: {len(movers2)}", flush=True)
    json.dump({'base': [str(base[0]), str(base[1])],
               'd1': {str(b): str(d1[b]) for b in control},
               'd2': {str(b): str(d2[b]) for b in control}},
              open('twist_deltas.json', 'w'))
    # linearity check on a few pairs
    mv = movers1[:6]
    import itertools
    ok = True
    with mp.Pool(6, initializer=init) as pool:
        pairs = list(itertools.combinations(mv, 2))[:8]
        for sb, (r1, r2) in pool.imap_unordered(worker, pairs):
            i, j = sb
            pred1 = base[0] + d1[i] + d1[j]
            pred2 = base[1] + d2[i] + d2[j]
            m = (pred1 == r1 and pred2 == r2)
            if not m: ok = False
            print(f"pair {sb}: linear={m}", flush=True)
    print(f"TWIST LINEAR IN BITS: {ok}", flush=True)

if __name__ == '__main__':
    main()
