#!/usr/bin/env python3
"""K35: do the 900 non-leaf free booleans matter?

Two things ride on this.  (1) If any of them is a real internal liveness bit, my drive()
seeds it to 0 and kills folds -- which would explain the B-half gap.  (2) If they are real
knobs at all, then Section 1's "the only boolean inputs are the 256 leaf selectors" is wrong
and the configuration space is larger than I claimed.

Test: drive the same ON-set three ways -- the 900 seeded 0, seeded 1, and left to the end so
the closure derives them if it can -- and compare the root halves against the composition."""
import sys, os, json
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from k26_drive import FORBID

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
handles, leafsel, otherbools, wires = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
S = FD.SHIFT
ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}


def run(on, mode):
    seed = {u: 0 for u in handles}
    for u in leafsel: seed[u] = 1 if u in on else 0
    if mode == 'zero':
        for u in otherbools: seed[u] = 0
        order = handles + leafsel + otherbools + defvars + wires
    elif mode == 'one':
        for u in otherbools: seed[u] = 1
        order = handles + leafsel + otherbools + defvars + wires
    else:                                    # derive-if-possible
        order = handles + leafsel + defvars + wires + otherbools
    v, _ = C.close(seed, order, forbid=FORBID)
    return (((v[12186] + S) % P, v[16742]), ((v[14853] + S) % P, v[24908]))


def f(es):
    R = FD.INF
    for e in es: R = FD.add(R, bypow[e])
    return R


for on_exp, ea, eb in [([3, 10], [], [3, 10]), ([0, 3, 10], [0], [3, 10]),
                       ([0, 1, 3], [0, 1], [3])]:
    on = set(exp2sel[e] for e in on_exp)
    print('\nON exponents', on_exp)
    for mode in ('zero', 'one', 'derive'):
        A, B = run(on, mode)
        am = (A == f(ea)) if ea else 'n/a'
        bm = (B == f(eb)) if eb else 'n/a'
        print('   otherbools=%-7s  A match %-5s  B match %-5s' % (mode, am, bm))
