#!/usr/bin/env python3
"""Propagate with a policy for branch points; report determined/free split."""
import pickle, os, collections, sys
from model import Model, load_assign
from prop import Engine
HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748

M = Model()
E = Engine(M)
policy = sys.argv[1] if len(sys.argv) > 1 else 'zero'
wit = None
if policy == 'wit':
    wit = load_assign(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))

val = [None] * NV
rounds = 0
while True:
    n, conf, br = E.propagate(val)
    rounds += 1
    known = sum(1 for x in val if x is not None)
    print(f"round {rounds}: +{n} known={known} conflicts={len(conf)} branch={len(br)}")
    if conf:
        cc = collections.Counter(k for k, _ in conf)
        print("   conflicts:", cc)
        for k, a in conf[:10]:
            print("     ", k, a, M.src[a][:110])
        break
    if not br:
        break
    # assign branch points
    changed = 0
    for u, a, roots in br:
        if val[u] is not None:
            continue
        if policy == 'wit':
            pick = wit[u] if wit[u] in roots else roots[0]
        else:
            pick = 0 if 0 in roots else roots[0]
        val[u] = pick
        changed += 1
    print(f"   assigned {changed} branch vars")
    if changed == 0:
        break

known = sum(1 for x in val if x is not None)
print(f"FINAL known {known}/{NV}, free {NV-known}")
pickle.dump(val, open(os.path.join(HERE, f'prop_{policy}.pkl'), 'wb'))

# score if we fill unknowns with 0
v2 = [0 if x is None else x for x in val]
fails, av, cv = M.eq_fail(v2)
print("score with unknowns=0:", M.ne - len(fails))
