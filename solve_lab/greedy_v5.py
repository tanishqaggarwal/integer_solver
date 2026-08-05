#!/usr/bin/env python3
"""Purpose-built greedy constraint-repair over bits, using the correct v5 model.

v5 forward-eval is a correct forward model (validated == 39,019 at all-0). For any
bit set it yields a deterministic assignment; count violated atoms. Greedily add
the bit that most reduces violations (candidates evaluated in parallel). Seed from
twist-fixing pairs. If it reaches 0, the witness is found; verify in Z."""
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

def fwd_modp(setbits):
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

def nviol_modp(val):
    A = _G['A']; n = 0
    for poly in A:
        s = 0
        for m, c in poly.items():
            t = c % P
            for x in m: t = (t*val[x]) % P
            s = (s+t) % P
        if s: n += 1
    return n

def eval_bits(setbits):
    return nviol_modp(fwd_modp(setbits))

def worker(args):
    cur, b = args
    return b, eval_bits(cur + [b])

def main():
    seed = [int(x) for x in sys.argv[1:]] or [710, 1858]
    init()
    control = json.load(open('control_bits.json'))
    cur = list(seed)
    base = eval_bits(cur)
    print(f"seed {cur}: {base} violations", flush=True)
    pool = mp.Pool(6, initializer=init)
    t0 = time.time()
    for rnd in range(40):
        cands = [(cur, b) for b in control if b not in cur]
        results = pool.map(worker, cands)
        best = min(results, key=lambda r: r[1])
        if best[1] >= base:
            print(f"round {rnd}: no improvement (min {best[1]} >= {base}); STUCK", flush=True)
            break
        cur.append(best[0]); base = best[1]
        print(f"round {rnd}: +bit {best[0]} -> {base} violations; set={sorted(cur)} ({time.time()-t0:.0f}s)", flush=True)
        if base == 0:
            init(); val = fwd_modp(cur)
            # exact Z verify
            from confluent_eval5 import make_forward
            solve = make_forward(_G['kind'], _G['info'], _G['seq'], _G['bestval'])
            vz = solve(list(_G['bestval']), cur)
            vio = sum(1 for poly in _G['A'] if any(True for _ in [0]) and 0)
            vio = 0
            for poly in _G['A']:
                s = 0
                for m, c in poly.items():
                    t = c
                    for x in m: t *= vz[x]
                    s += t
                if s: vio += 1
            print(f"  Z-verify: {vio} violations", flush=True)
            if vio == 0:
                json.dump({f"x_{i}": vz[i] for i in range(NVARS)}, open('cand_SOLVED.json', 'w'))
                print("  *** SOLVED! wrote cand_SOLVED.json ***", flush=True)
            break
    pool.close()

if __name__ == '__main__':
    main()
