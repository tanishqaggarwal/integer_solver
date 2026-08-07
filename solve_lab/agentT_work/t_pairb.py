#!/usr/bin/env python3
"""AUDIT T31b -- is the two-wire |S|=8 closure true OF THE INSTANCE, not just of L's 9,032-atom
engine?  Same test as t_S2b: evaluate all 39,033 atoms of F's certified-faithful decomposition at
the dumped assignment and check the nonzero set is exactly the two target congruences and that
their equation footprint is exactly checker.py's failing set.
(F's circ4.pkl was wiped with every other *.pkl by the restart; rebuilt in agentT_work/mirror/F.)"""
import os, sys, json, pickle, collections
T = os.path.dirname(os.path.abspath(__file__)); LAB = os.path.abspath(os.path.join(T, '..'))
F = os.path.join(T, 'mirror', 'F'); sys.path.insert(0, F); sys.path.insert(0, LAB)
from fwd import compile_node
import checker as CK
d = pickle.load(open(os.path.join(F, 'circ4.pkl'), 'rb'))
atoms = d['atoms']; eqrows = d['eqrows']; names = list(atoms); idx = {a: i for i, a in enumerate(names)}
a2e = collections.defaultdict(set)
for e, row in enumerate(eqrows):
    for k, a in row:
        a2e[idx[a]].add(e)
prog = compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']', '<at>', 'exec')
NV = 38748
codes, varsets = CK.load_equations()
print('F parse: %d atoms, %d equations' % (len(names), len(eqrows)))
targets = sys.argv[1:] or [os.path.join(T, 'close_T8pair.json')]
for fn in targets:
    v = [0]*NV
    for k, val in json.load(open(fn)).items():
        v[int(k[2:])] = int(val)
    r = [0]*len(names); exec(prog, {'v': v, 'r': r, '__builtins__': {}})
    nz = [i for i in range(len(names)) if r[i]]
    fails = CK.evaluate_all(codes, v)
    foot = set()
    for i in nz:
        foot |= a2e[i]
    print('\n== %s ==' % os.path.basename(fn))
    print('   nonzero atoms in F\'s parse : %d' % len(nz))
    for i in nz:
        print('        %s' % names[i][:90])
    print('   checker satisfied          : %d/%d  (%d failing)' % (len(eqrows)-len(fails),
                                                                   len(eqrows), len(fails)))
    print('   failing                    : %s' % fails)
    print('   equation footprint         : %d' % len(foot))
    print('   footprint == failing set ? %s' % (foot == set(fails)))
    if foot != set(fails):
        print('     in footprint not failing (cancellation): %s' % sorted(foot-set(fails)))
        print('     failing not in footprint (UNEXPLAINED) : %s' % sorted(set(fails)-foot))
