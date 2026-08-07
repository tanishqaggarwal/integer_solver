#!/usr/bin/env python3
"""INVERTED SEARCH: rank footprints by how cheap they are, BEFORE asking reachability.

Every search in this campaign went configuration-first and priced the footprint it landed in.
This goes the other way: enumerate candidate live-atom supports S, price each by
    cost(S) = |equations touched by S|  -  max #{touched equations a nonzero value vector kills}
and rank.  cost(S) is the failing floor for any assignment whose nonzero-atom set is exactly S.
Beating 39,026 needs a support with cost(S) <= 6.
"""
import sys, json, collections, itertools
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import E

NA = len(E.res)
atom_eqs = collections.defaultdict(list)     # atom -> [(eq, coeff)]
for j, rr in enumerate(E.eqres):
    for c, i in rr: atom_eqs[i].append((j, c))
print('atoms=%d  equations=%d  incidences=%d' % (NA, len(E.eqres), sum(len(v) for v in atom_eqs.values())), flush=True)

# ---- |S| = 1 : cost = number of equations containing the atom with nonzero coeff ----
occ = {i: len({j for j, c in v if c}) for i, v in atom_eqs.items()}
hist = collections.Counter(occ.values())
print('single-atom footprint cost distribution (cost -> how many atoms):')
for k in sorted(hist)[:20]: print('   cost %-4d : %d atoms' % (k, hist[k]), flush=True)
best1 = sorted(occ.items(), key=lambda kv: kv[1])[:40]
print('cheapest single atoms:', best1[:15], flush=True)
json.dump({'occ_hist': {str(k): v for k, v in sorted(hist.items())},
           'cheapest_single': [[int(a), int(c)] for a, c in best1]},
          open('runs/footprints1.json', 'w'), indent=1)
