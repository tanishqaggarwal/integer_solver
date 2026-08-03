#!/usr/bin/env python3
"""Genetic algorithm over the 255 control bits, minimizing mod-P forward-eval
violation count. Population + tournament selection + uniform crossover + mutation,
biased to activate x_12779 (the slack gate). Long-running; checkpoints improvements."""
import json, time
from confluent_eval5 import build5, make_forward
from propagate import NVARS
from modp import P

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
    control = json.load(open('control_bits.json')); nc = len(control)

    def nviol(mask):
        sb = [control[i] for i in range(nc) if mask[i]]
        val = solveP(list(bm), sb); n = 0
        for poly in A:
            s = 0
            for m, c in poly.items():
                t = c % P
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            if s: n += 1
        return n

    st = 20240711
    def rnd():
        nonlocal st; st = (st*6364136223846793005+1442695040888963407) & ((1<<64)-1); return st >> 33
    POP = 24
    def randmask(k):
        m = [0]*nc
        for _ in range(k):
            m[rnd() % nc] = 1
        return m
    pop = [randmask(1+rnd() % 8) for _ in range(POP)]
    fit = [nviol(m) for m in pop]
    best = min(fit); print(f"init best {best} ({time.time()-t0:.0f}s)", flush=True)
    gen = 0
    while time.time()-t0 < 9000:  # up to 2.5h
        gen += 1
        # tournament select 2 parents
        def sel():
            i, j = rnd() % POP, rnd() % POP
            return pop[i] if fit[i] <= fit[j] else pop[j]
        p1, p2 = sel(), sel()
        child = [p1[i] if rnd() & 1 else p2[i] for i in range(nc)]
        # mutation
        for _ in range(1 + rnd() % 3):
            child[rnd() % nc] ^= 1
        fc = nviol(child)
        # replace worst
        wi = max(range(POP), key=lambda i: fit[i])
        if fc < fit[wi]:
            pop[wi] = child; fit[wi] = fc
        if fc < best:
            best = fc
            sb = sorted(control[i] for i in range(nc) if child[i])
            print(f"  gen {gen}: NEW BEST {best} bits={sb[:20]} ({time.time()-t0:.0f}s)", flush=True)
            if best < 4:
                # Z verify
                solveZ = make_forward(kind, info, seq, bestval)
                vz = solveZ(list(bestval), sb); vio = 0
                for poly in A:
                    s = 0
                    for m, c in poly.items():
                        t = c
                        for x in m: t *= vz[x]
                        s += t
                    if s: vio += 1
                print(f"    Z-verify: {vio}", flush=True)
                if vio == 0:
                    json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json', 'w'))
                    print("    *** SOLVED ***", flush=True); return
                else:
                    json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open(f'ga_best_{vio}.json', 'w'))
        if gen % 500 == 0:
            print(f"  gen {gen}: best {best}, pop-min {min(fit)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"GA done: best {best} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
