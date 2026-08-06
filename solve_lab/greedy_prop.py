#!/usr/bin/env python3
"""Custom greedy descent over bits via propagation (no SAT/SMT).

Start from a twist-fixing pair (zeros the 3 twist checks but cascades to ~20
violated atoms). Greedily add the control bit that most reduces total violated
atoms (each candidate evaluated by full direction-agnostic propagation). If the
cascade shrinks toward 0 the witness is being assembled bit by bit."""
import json, time, sys
from propagate import load_atoms
from solve_forward import Engine
NVARS = 38748

atoms = load_atoms()
control = json.load(open('control_bits.json'))
BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def viol_count(val):
    n = 0
    for poly in atoms:
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: n += 1
    return n

def prop_solve(setbits):
    eng = Engine(atoms); eng.propagate()
    for b in setbits:
        if eng.val[b] is None: eng.assign(b, 1); eng.propagate()
    for v in [x for x in control if eng.val[x] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    val = [x if x is not None else 0 for x in eng.val]
    return viol_count(val), val

def main():
    t0 = time.time()
    cur = set([710, 1858])
    base, _ = prop_solve(list(cur))
    print(f"start {sorted(cur)}: {base} violations ({time.time()-t0:.0f}s)", flush=True)
    # candidate pool: 232-movers + 22 bits (all control), avoid current
    for rnd in range(20):
        best = (base, None)
        for b in control:
            if b in cur: continue
            v, _ = prop_solve(list(cur | {b}))
            if v < best[0]: best = (v, b)
        if best[1] is None:
            print(f"round {rnd}: no improving bit (min stays {base}) ({time.time()-t0:.0f}s)", flush=True)
            break
        cur.add(best[1]); base = best[0]
        print(f"round {rnd}: added bit {best[1]} -> {base} violations, set={sorted(cur)} ({time.time()-t0:.0f}s)", flush=True)
        if base == 0:
            v, val = prop_solve(list(cur))
            json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_greedy.json', 'w'))
            print("*** 0 violations! wrote cand_greedy.json ***", flush=True)
            break

if __name__ == '__main__':
    main()
