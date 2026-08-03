#!/usr/bin/env python3
"""Empirically map the twist quadruple (x_9770,x_18274,x_3183,x_17728) under the
evaluator over many bit settings. Central question: can x_3183==x_17728 with both
nonzero (the one HARD twist constraint)? Also record x_20510*x_31302 vs
6033033*(x_18274-x_9770) (the a1817 slack balance) and x_9897*x_12518 (a1818,
must be 0). Separates 22-side bits from 233-side bits to see each side's image."""
import json, time
from confluent_eval5 import build5, make_forward

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    b22 = [b for b in control if b in BITS22]
    b233 = [b for b in control if b not in BITS22]
    print(f"22-side bits: {len(b22)}, 233-side bits: {len(b233)}", flush=True)

    st = 12345
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    def rec(bits):
        v = solve(list(bestval), bits)
        return (v[9770], v[18274], v[3183], v[17728], v[20510]*v[31302], 6033033*(v[18274]-v[9770]), v[9897]*v[12518])

    base = rec([])
    print(f"all-0: 9770={base[0]} 18274={base[1]} 3183={base[2]} 17728={base[3]}", flush=True)
    print(f"       x20510*x31302={base[4]} 6033033*(18274-9770)={base[5]} x9897*x12518={base[6]}", flush=True)

    # image of x_3183 over 22-side; image of x_17728 over 233-side
    img3183 = set(); img17728 = set(); img9770 = set(); img18274 = set()
    for _ in range(400):
        k = 1 + rnd() % 6
        s22 = [b22[rnd() % len(b22)] for _ in range(k)]
        v = solve(list(bestval), s22)
        img3183.add(v[3183]); img9770.add(v[9770])
    for _ in range(600):
        k = 1 + rnd() % 8
        s233 = [b233[rnd() % len(b233)] for _ in range(k)]
        v = solve(list(bestval), s233)
        img17728.add(v[17728]); img18274.add(v[18274])
    print(f"\n|img x_3183 (22-side)|={len(img3183)}  |img x_9770 (22-side)|={len(img9770)}", flush=True)
    print(f"|img x_17728 (233-side)|={len(img17728)}  |img x_18274 (233-side)|={len(img18274)}", flush=True)
    inter = img3183 & img17728
    print(f"x_3183(22) ∩ x_17728(233) = {len(inter)} common values: {sorted(inter)[:6]}", flush=True)
    # does mixing sides let x_3183 move with 233-bits or x_17728 with 22-bits?
    mv3183_by233 = set(); mv17728_by22 = set()
    for b in b233:
        mv3183_by233.add(solve(list(bestval), [b])[3183])
    for b in b22:
        mv17728_by22.add(solve(list(bestval), [b])[17728])
    print(f"x_3183 distinct under single 233-bit: {len(mv3183_by233)} (base {base[2]})", flush=True)
    print(f"x_17728 distinct under single 22-bit: {len(mv17728_by22)} (base {base[3]})", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
