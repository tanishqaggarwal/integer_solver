#!/usr/bin/env python3
"""Re-score every assignment I generated with the LAB'S checker.py, not agent F's
evaluator E.  E was measured to over-report failures on the deliverable (13 vs 7)."""
import sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
from cfgscan import run_cfg, E
from fwd import NV
import checker

codes, _ = checker.load_equations()
print('equations', len(codes), flush=True)
out = {}
def score(v, tag):
    f = checker.evaluate_all(codes, list(v))
    r = E.run(list(v)); fe = sorted(E.score(r))
    out[tag] = {'checker_failing': len(f), 'checker_score': len(codes) - len(f),
                'E_failing': len(fe), 'E_score': len(codes) - len(fe),
                'E_overreport': len(fe) - len(f),
                'first_failing': sorted(f)[:12]}
    print('%-22s checker %d failing -> score %d   |  E %d failing -> %d   (E over by %d)'
          % (tag, len(f), len(codes) - len(f), len(fe), len(codes) - len(fe), len(fe) - len(f)), flush=True)
    json.dump(out, open('runs/rescore.json', 'w'), indent=1)
    if len(codes) - len(f) > 39026:
        json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]},
                  open('runs/BEAT_%d.json' % (len(codes) - len(f)), 'w'))
        print('*** BEATS BASELINE ***', flush=True)
    return len(codes) - len(f)

d = checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
score(d, 'deliverable')
for cfg in ([24601], [2081], [47], [24601, 2081], []):
    try:
        sc, nz, ok, v = run_cfg(cfg)
        score(v, 'gs2 ' + (','.join(map(str, cfg)) or '<empty>'))
    except Exception as e:
        print('gs2', cfg, 'ERR', repr(e)[:120], flush=True)
