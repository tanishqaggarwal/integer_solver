#!/usr/bin/env python3
"""agent AF, step 10: recover the 383 law blocks and their liveness gates from the conditions."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import (atoms, defs, defc, val, lift, Pval, find, pp, expand, varsof, shape_of)
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
conds = C['conds']; sel = T['sel']; selset = set(sel)

def strip_k(n):
    """peel integer factors: return (k, node)"""
    k = 1
    while True:
        if n[0] == 'neg':
            k = -k; n = n[1]; continue
        if n[0] == '*':
            c1 = is_const(n[1]); c2 = is_const(n[2])
            if c1 is not None:
                k *= c1; n = n[2]; continue
            if c2 is not None:
                k *= c2; n = n[1]; continue
        return k, n

prod = []      # (cond_index, k, A, B)
diff = []      # (cond_index, k, X, Y, sign)
pins = []
other = []
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    k, n = strip_k(Ex)
    if n[0] == '*' and n[1][0] == 'v' and n[2][0] == 'v':
        prod.append((i, k, find(n[1][1]), find(n[2][1])))
    elif n[0] == '*' and n[1][0] == 'v' and n[2][0] == '-' and n[2][1][0] == 'v' \
         and is_const(n[2][2]) is not None:
        pins.append((i, k, find(n[1][1]), find(n[2][1][1]), is_const(n[2][2])))
    elif n[0] in ('+', '-') and n[1][0] == 'v' and n[2][0] == 'v':
        diff.append((i, k, find(n[1][1]), find(n[2][1]), 1 if n[0] == '-' else -1))
    else:
        other.append((i, k, n))
print('condition Expr shapes:  product %d | pin %d | difference %d | other %d'
      % (len(prod), len(pins), len(diff), len(other)))
print('other shapes:', Counter(shape_of(n) for (_, _, n) in other).most_common(6))

fac = Counter()
for (i, k, A, B) in prod:
    fac[A] += 1; fac[B] += 1
rep = [v for v, n in fac.items() if n >= 2]
print('variables appearing as a factor in >=2 product conditions: %d' % len(rep))
print('  multiplicity histogram:', dict(Counter(fac[v] for v in rep)))

# a gate is a factor whose complement (1-g) also appears as a factor
notof = {}
for r, dl in defc.items():
    for aid, rhs in dl:
        if rhs[0] == '-' and is_const(rhs[1]) == 1 and rhs[2][0] == 'v':
            notof.setdefault(r, find(rhs[2][1]))
gates = {}
for v in fac:
    if v in notof and notof[v] in fac:
        gates[notof[v]] = v      # gates[L] = var holding (1-L)
print('gate/complement pairs both used as condition factors: %d' % len(gates))
mult = Counter()
for L, nL in gates.items():
    mult[(fac[L], fac[nL])] += 1
print('  (uses of L, uses of 1-L) histogram:', dict(mult))

pickle.dump({'prod': prod, 'pins': pins, 'diff': diff, 'other': other,
             'gates': gates, 'notof': notof, 'fac': dict(fac)},
            open(os.path.join(HERE, 'af_blocks.pkl'), 'wb'), 2)
