#!/usr/bin/env python3
"""Long-running simulated-annealing campaign over the 255 control bits, minimizing
the TRUE forward-eval violation count (mod-P proxy; Z-verify any improvement).
Multi-restart, pair-flips, deterministic pseudo-random (no RNG: LCG on an int seed
passed via steps). Checkpoints any state with < 4 violated atoms. all-0 is the known
local min at 4; this hunts for the isolated 0-violation witness."""
import json, time, sys
from confluent_eval5 import build5, make_forward
from propagate import NVARS
from modp import P, inv

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    control = json.load(open('control_bits.json'))
    nc = len(control)

    Aatoms = A
    def nviol(setbits):
        val = solveP(list(bm), list(setbits))
        n = 0
        for poly in Aatoms:
            s = 0
            for m, c in poly.items():
                t = c % P
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            if s: n += 1
        return n

    # LCG pseudo-random
    st = 12345
    def rnd():
        nonlocal st
        st = (st*1103515245 + 12345) & 0x7fffffff
        return st

    best_global = 4
    n_eval = 0
    n_restart = 0
    T0 = 20.0
    import math
    while time.time() - t0 < 5400:   # up to 90 min
        n_restart += 1
        # start: mostly from all-0 with a few random bits, occasionally denser
        k = 1 + (rnd() % (3 if n_restart % 2 else 12))
        cur = set()
        for _ in range(k):
            cur.add(control[rnd() % nc])
        cur_v = nviol(sorted(cur)); n_eval += 1
        T = T0
        for step in range(300):
            T *= 0.985
            b = control[rnd() % nc]
            new = set(cur)
            if b in new: new.discard(b)
            else: new.add(b)
            # occasional pair flip
            if rnd() % 3 == 0:
                b2 = control[rnd() % nc]
                if b2 in new: new.discard(b2)
                else: new.add(b2)
            nv = nviol(sorted(new)); n_eval += 1
            d = nv - cur_v
            if d <= 0 or (T > 0.01 and rnd()/0x7fffffff < math.exp(-d/T)):
                cur, cur_v = new, nv
                if cur_v < best_global:
                    best_global = cur_v
                    print(f"  NEW BEST mod-P viol {cur_v} at {sorted(cur)} ({time.time()-t0:.0f}s)", flush=True)
                    # Z-verify
                    solveZ = make_forward(kind, info, seq, bestval)
                    vz = solveZ(list(bestval), sorted(cur))
                    vio = sum(1 for poly in Aatoms if any(True for _ in [0]) and 0)
                    vio = 0
                    for poly in Aatoms:
                        s = 0
                        for m, c in poly.items():
                            t = c
                            for x in m: t *= vz[x]
                            s += t
                        if s: vio += 1
                    print(f"    Z-verify: {vio} violated", flush=True)
                    if vio < 4:
                        json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open(f'sa_best_{vio}.json','w'))
                        print(f"    *** improved partial saved (Z viol {vio}) ***", flush=True)
                    if vio == 0:
                        json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json','w'))
                        print("    *** SOLVED! ***", flush=True); return
        if n_restart % 20 == 0:
            print(f"  ...{n_restart} restarts, {n_eval} evals, best {best_global} ({time.time()-t0:.0f}s)", flush=True)
    print(f"campaign done: {n_restart} restarts, {n_eval} evals, best mod-P viol {best_global} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
