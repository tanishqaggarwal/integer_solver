#!/usr/bin/env python3
"""K37: test the PREMISE my partition theorem rests on.

The theorem needs: at each stage, the two slot values are the compositions of the live leaves
in their supports.  Two ways that can fail, and Q/T's "routing is a constraint, not a
propagation" bears on both:

  TEST 1 (identification).  Maybe the composition IS on a wire and I named the wrong wire.
  For ON={e3,e10} scan EVERY variable for the predicted composition's coordinates.

  TEST 2 (forcedness).  Maybe the constraints do not force the slot at all.  Seed the slot
  wire to a random wrong value BEFORE closing, so the closure cannot derive it, and see
  whether the system still closes with zero nonzero atoms.  Closes => the wire is free and
  routing does not pin it => my premise is false as stated.  Nonzero atoms => it is pinned.
"""
import sys, os, json, random
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
ORDER = handles + leafsel + otherbools + defvars + wires
S = FD.SHIFT
ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}


def close(on, pre=None):
    seed = {u: 0 for u in handles}
    for u in leafsel: seed[u] = 1 if u in on else 0
    for u in otherbools: seed[u] = 0
    order = ORDER
    if pre:
        for u, val in pre.items(): seed[u] = val
        order = handles + list(pre) + leafsel + otherbools + defvars + wires
    v, _ = C.close(seed, order, forbid=FORBID)
    return v


def comp(es):
    R = FD.INF
    for e in es: R = FD.add(R, bypow[e])
    return R


print('=' * 70)
print('TEST 1 - is the predicted composition anywhere on any wire?')
for ON_EXP in ([3, 10], [3, 5], [0, 1]):
    on = set(exp2sel[e] for e in ON_EXP)
    v = close(on)
    pr = comp(ON_EXP)
    hx = [u for u in range(NV) if (v[u] + S) % P == pr[0]]
    hxr = [u for u in range(NV) if v[u] == pr[0]]
    hy = [u for u in range(NV) if v[u] == pr[1]]
    print('ON exps %-8s composition = (%d..., %d...)' % (ON_EXP, pr[0] % 10 ** 9, pr[1] % 10 ** 9))
    print('   wires holding X-shift : %s' % (hx[:8] or 'NONE'))
    print('   wires holding X raw   : %s' % (hxr[:8] or 'NONE'))
    print('   wires holding Y       : %s' % (hy[:8] or 'NONE'))

print()
print('=' * 70)
print('TEST 2 - are the root slot wires FORCED, or free?')
random.seed(11)
on = set(exp2sel[e] for e in (0, 1, 3))
base = close(on)
nz0 = len(C.nzatoms(base))
print('baseline ON={e0,e1,e3}: nonzero atoms mod p =', nz0)
for w, nm in ((12186, 'root A.x'), (16742, 'root A.y'), (14853, 'root B.x'), (24908, 'root B.y')):
    bad = random.randrange(2, P - 2)
    v = close(on, pre={w: bad})
    nz = len(C.nzatoms(v))
    kept = (v[w] == bad)
    print('  x%-6d %-9s seeded to a random value: kept=%s  nonzero atoms=%d  %s'
          % (w, nm, kept, nz, 'FREE (routing does NOT pin it)' if kept and nz <= nz0
             else 'pinned / rejected'))
