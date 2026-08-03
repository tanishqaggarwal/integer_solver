#!/usr/bin/env python3
"""Test linearity of x_9770,x_3183 in the 22 bits and x_18274,x_17728,x_26977,
x_9982 in the 233 bits. Compute single-flip deltas; verify additivity on random
subsets. If linear, the whole twist is a linear system over binary vars."""
import json, time
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]
    base = solve(list(bestval), [])
    WATCH = [9770, 3183, 18274, 17728, 26977, 9982]
    b0 = {w: base[w] for w in WATCH}
    print("base:", {f'x_{w}': b0[w] for w in WATCH}, flush=True)

    # single-flip deltas
    delta = {w: {} for w in WATCH}
    for b in control:
        val = solve(list(bestval), [b])
        for w in WATCH:
            dv = val[w] - b0[w]
            if dv: delta[w][b] = dv
    for w in WATCH:
        m22 = [b for b in delta[w] if b in set(BITS22)]
        m233 = [b for b in delta[w] if b not in set(BITS22)]
        print(f"x_{w}: moved by {len(m22)} of 22, {len(m233)} of 233", flush=True)

    # linearity: test additivity on subsets
    def test(w, bitpool, sizes):
        ok = True
        for sz in sizes:
            S = bitpool[:sz]
            val = solve(list(bestval), S)
            actual = val[w] - b0[w]
            pred = sum(delta[w].get(b, 0) for b in S)
            lin = (actual == pred)
            print(f"  x_{w} |S|={sz}: {'LIN' if lin else 'NONLIN diff='+str(actual-pred)}", flush=True)
            if not lin: ok = False
        return ok
    print("\n-- x_9770 over 22 bits --")
    p22 = [b for b in BITS22 if b in delta[9770]] + [b for b in BITS22 if b not in delta[9770]]
    test(9770, p22, [2, 5, 10, 22])
    print("-- x_3183 over 22 bits --")
    p22b = [b for b in BITS22 if b in delta[3183]] + [b for b in BITS22 if b not in delta[3183]]
    test(3183, p22b, [2, 5, 10, 22])
    print("-- x_18274 over 233 bits --")
    p233 = [b for b in bits233 if b in delta[18274]] + [b for b in bits233 if b not in delta[18274]]
    test(18274, p233, [2, 5, 20, 60, 120, 211])
    print("-- x_17728 over 233 bits --")
    p233b = [b for b in bits233 if b in delta[17728]] + [b for b in bits233 if b not in delta[17728]]
    test(17728, p233b, [2, 5, 20, 60, 120, 211])
    print("-- x_26977 (target of atom1817 3rd var) over ALL --")
    pall = list(delta[26977].keys()) + [b for b in control if b not in delta[26977]]
    test(26977, pall, [2, 5, 20])
    print("-- x_9982 over ALL --")
    pall2 = list(delta[9982].keys()) + [b for b in control if b not in delta[9982]]
    test(9982, pall2, [2, 5, 20])

    # overlap of x_18274 and x_17728 bit supports
    s1 = set(delta[18274]); s2 = set(delta[17728])
    print(f"\nx_18274 support {len(s1)}, x_17728 support {len(s2)}, overlap {len(s1&s2)}", flush=True)

    json.dump({'base': {str(w): str(b0[w]) for w in WATCH},
               'delta': {str(w): {str(b): str(v) for b, v in delta[w].items()} for w in WATCH}},
              open('watch_deltas.json', 'w'))
    print(f"wrote watch_deltas.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
