"""Step 12: equation cost of each residual atom at stateA2 / stateA3."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
d = pickle.load(open('atoms.pkl', 'rb'))
eq_terms = d['eq_terms']

for st in ('quad/stateA2.json', 'quad/stateA3.json'):
    v = H.load_assignment(st)
    live = sorted(a for a in range(len(polys)) if evalpoly(polys[a], v) != 0)
    av = {a: evalpoly(polys[a], v) for a in live}
    cand = eqs_of(live)
    ff = H.evaluate(CODES, v, cand)
    print(f'\n=== {st}: live atoms {live}  FAIL={len(ff)}')
    # per failing equation: which live atoms it contains and with what coeffs
    cnt = collections.Counter()
    for i in ff:
        m, sq, tl = eq_terms[i]
        cc = {a: c for c, a in tl if a in av}
        cnt[frozenset(cc)] += 1
    for k, n in cnt.most_common():
        print(f'   {n:3d} eqs involve atoms {sorted(k)}')
    # cost of each single atom alone
    for a in live:
        eqs = set(atom2eq.get(a, []))
        solo = [i for i in ff if set(x for c, x in eq_terms[i][2]) & set(live) <= {a}]
        print(f'   atom {a}: appears in {len(eqs)} eqs; failing eqs where it is the ONLY live atom: {len(solo)}')
