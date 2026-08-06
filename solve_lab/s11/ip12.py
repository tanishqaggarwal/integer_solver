"""IP #12 -- is the p-factor of the invariant UNIVERSAL?

At the checkpoint the invariant is 2458959*p; at the session's 39,018 state it is 8640431*p.
Compute the invariant D at every saved state and test whether p | D always.  If so, that is a
clean universal statement of the trapdoor in integer-programming language:

    the residual integer program is always consistent over Q, and its invariant factor
    always contains p.
"""
import sys, os, json, math, time, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip8 import build
from ip7 import load_raw
from ip9 import rational_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)

STATES = ['../best/new_instance_partial_39026.json',
          'data/finish3_named.json', 'data/closehit2.json', 'data/three.json',
          'data/quad3_hit.json', 'data/tri7_best.json', 'data/uv01_full.json',
          'data/eqopt2_named.json', 'data/cheapdefect_named.json']

rows = []
for rel in STATES:
    path = os.path.join(HERE, rel)
    if not os.path.exists(path):
        continue
    try:
        v = load_raw(path)
        v, FAIL, used, M, rhs, nf = build(v, verbose=False)
        if not M:
            continue
        t0 = time.time()
        x, piv = rational_solve(M, rhs)
        if x is None:
            print(f"{os.path.basename(rel):32s} failing={nf:3d}  INCONSISTENT over Q")
            rows.append((rel, nf, None, None))
            continue
        D = 1
        for t in x:
            D = D * t.denominator // math.gcd(D, t.denominator)
        cof = D // P if D % P == 0 else None
        print(f"{os.path.basename(rel):32s} failing={nf:3d}  D={len(str(D)):3d} digits  "
              f"p | D : {D % P == 0}   D/p = {cof if cof and cof < 10**12 else (str(cof)[:18]+'..' if cof else '-')}"
              f"   ({time.time()-t0:.0f}s)", flush=True)
        rows.append((rel, nf, D, cof))
    except Exception as e:
        print(f"{os.path.basename(rel):32s} error {e}")

ok = [r for r in rows if r[2] is not None]
print()
print(f"states analysed: {len(ok)}")
print(f"  consistent over Q          : {len(ok)}/{len(ok)}")
print(f"  invariant divisible by p   : {sum(1 for r in ok if r[2] % P == 0)}/{len(ok)}")
print(f"  cofactors D/p              : {sorted({r[3] for r in ok if r[3]})}")
