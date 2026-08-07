#!/usr/bin/env python3
"""Try to REALIZE the 39,027 floor: relax selectors x33095 + x19326 (union of their
boolean-ness equations = 6) and let gs2's repair choose their values.
Nothing is claimed unless checker.py verifies it."""
import sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import gs2
from cfgscan import PIN, PARTNER, TREE
from fwd import NV
import checker

codes, _ = checker.load_equations()
print('equations', len(codes), flush=True)
best = 39026
def run(tag, bits, relax, vals):
    v = [0] * NV
    fr = set(PIN)
    for k, x in PIN.items(): v[k] = x
    for b in bits:
        v[b] = 1; fr.add(b)
        pa = PARTNER[TREE[b]]; v[pa] = 1; fr.add(pa)
    for r, val in zip(relax, vals):
        v[r] = val; fr.add(r)          # frozen at a NON-BOOLEAN value
    try:
        v2, ok = gs2.solve(v, verbose=False, frozen=set(fr))
    except Exception as e:
        print('%-40s gs2 ERR %s' % (tag, repr(e)[:80]), flush=True); return
    f = checker.evaluate_all(codes, list(v2))
    sc = len(codes) - len(f)
    print('%-40s checker score %d  (%d failing)  ok=%s' % (tag, sc, len(f), ok), flush=True)
    global best
    if sc > best:
        best = sc
        json.dump({'x_%d' % i: v2[i] for i in range(NV) if v2[i]}, open('runs/BEAT_%d.json' % sc, 'w'))
        print('*** BEATS BASELINE, written to runs/BEAT_%d.json ***' % sc, flush=True)
    return sc

P = 2 ** 256  # placeholder magnitude for a non-boolean value
for vals in ([2, 2], [2, 3], [7376877, 5113045]):
    run('relax x33095,x19326 = %s' % vals, [24601, 2081], [33095, 19326], vals)
for vals in ([2], [3]):
    run('relax x33095 = %s' % vals, [24601, 2081], [33095], vals)
    run('relax x33095 = %s (cfg 24601)' % vals, [24601], [33095], vals)
run('baseline cfg 24601,2081 (no relax)', [24601, 2081], [], [])
