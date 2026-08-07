#!/usr/bin/env python3
"""Do agent F's evaluator E and the lab's checker.py agree on the SAME assignment?
Everything I reported was on E's scale; this prices that assumption."""
import sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
from cfgscan import E
from fwd import NV
import checker

v = checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
t = time.time(); codes, varsets = checker.load_equations(); print('eqs loaded', len(codes), '%.0fs' % (time.time() - t))
t = time.time(); f_chk = checker.evaluate_all(codes, v); print('checker failing', len(f_chk), sorted(f_chk)[:12], '%.0fs' % (time.time() - t))
t = time.time(); r = E.run(list(v)); f_E = sorted(E.score(r)); print('E       failing', len(f_E), f_E[:12], '%.0fs' % (time.time() - t))
print('E-only :', sorted(set(f_E) - set(f_chk))[:20])
print('chk-only:', sorted(set(f_chk) - set(f_E))[:20])
json.dump({'checker': sorted(f_chk), 'E': f_E}, open('runs/crosscheck.json', 'w'))
