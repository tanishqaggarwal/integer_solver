"""FAST p-adic solvability test by Hensel lifting.

M x = rhs is solvable over Z_p iff it lifts to every power of p.  Writing x = x0 + p*x1 + ...:

    step 0 :  M x0 = rhs           (mod p)          -- a GF(p) solve
    step k :  M xk = (rhs - M*(partial))/p^k  (mod p) -- another GF(p) solve

Each step is fast, so this decides the p-part of the invariant in seconds instead of the hours
an integer HNF takes -- which is what has been throttling every experiment.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip14 import gf_solve
from ip7 import load_raw
from closed import build_closed
from ip8 import build as build_small
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def hensel(M, rhs, q=P, steps=6, verbose=True):
    """returns the number of successful lifting steps (steps = fully p-adically solvable so far)"""
    m = len(M)
    n = len(M[0])
    acc = [0] * n
    cur = list(rhs)
    for k in range(steps):
        x = gf_solve(M, [c % q for c in cur], q)
        if x is None:
            if verbose:
                print(f"     lift fails at p^{k+1}")
            return k
        for j in range(n):
            acc[j] += x[j] * (q ** k)
        nxt = []
        for i in range(m):
            r = cur[i] - sum(M[i][j] * x[j] for j in range(n))
            if r % q:
                if verbose:
                    print(f"     residual not divisible by p at step {k}")
                return k
            nxt.append(r // q)
        cur = nxt
        if all(c == 0 for c in cur):
            if verbose:
                print(f"     EXACT integer solution found after {k+1} steps")
            return steps + 1
    return steps


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    print(f"=== {os.path.basename(src)}")
    for cap in [130, 300, 500, 900, 1400]:
        v = load_raw(src)
        t0 = time.time()
        if cap == 130:
            v2, FAIL, used, M, rhs, nf = build_small(v, verbose=False)
        else:
            v2, FAIL, ROWS, used, M, rhs = build_closed(v, maxrows=cap, verbose=False)
        k = hensel(M, rhs, verbose=False)
        tag = ('SOLVABLE over Z' if k > 6 else f'lifts {k} step(s) then fails')
        print(f"  cap {cap:5d}: system {len(M)}x{len(M[0])}  ->  {tag}   ({time.time()-t0:.0f}s)",
              flush=True)
