"""Test EVERY single-row drop (failing AND collateral).

ip8 only ever dropped rows from the FAILING block.  But a feasible subsystem may need to
sacrifice a currently-SATISFIED equation instead -- and that is exactly what the checkpoint
itself does (it breaks five gates).  If dropping one row makes the whole 130-row system
integer-solvable, the result is ONE failing equation: score 39,032.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip7 import atomval, load_raw
from ip8 import build
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)

LAB = os.path.join(HERE, '..')
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v0 = load_raw(src)
print("===", os.path.basename(src), flush=True)
v, FAIL, used, M, rhs, nf = build(v0)
m = len(M)
print(f"  testing all {m} single-row drops (rows 0..{nf-1} are the failing ones)", flush=True)
t0 = time.time()
hits = []
for i in range(m):
    keep = [j for j in range(m) if j != i]
    x = solve_int([M[j] for j in keep], [rhs[j] for j in keep])
    if x is not None:
        hits.append((i, x))
        print(f"  *** ROW {i} DROP WORKS ({'FAILING' if i < nf else 'collateral'}) "
              f"({time.time()-t0:.0f}s)", flush=True)
        # apply and measure
        vv = [t for t in v]
        for j, u in enumerate(used):
            vv[u] += x[j]
        AV = [atomval(a, vv) for a in range(L.NA)]
        F = [e for e in range(L.NEQ)
             if sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) != 0]
        print(f"      applied -> failing={len(F)} score={L.NEQ-len(F)}", flush=True)
        if len(F) < 7:
            json.dump({('x_%d' % t): vv[t] for t in range(L.NVARS)},
                      open(os.path.join(HERE, 'data', f'drop1_row{i}.json'), 'w'))
            print("      SAVED", flush=True)
        break
    if i % 10 == 0:
        print(f"    {i}/{m} ({time.time()-t0:.0f}s)", flush=True)
if not hits:
    print(f"  no single-row drop works ({time.time()-t0:.0f}s)")
