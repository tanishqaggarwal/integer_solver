#!/usr/bin/env python3
"""KILL-TEST for my accumulator model.

My solve2.py assumes a chain with an accumulator seeded at L_0 whose intermediate values are
group elements carried on wires, gated by selector variables.  T and Q measured that routing is a
CONSTRAINT not a propagation, and that the deliverable keeps live leaves with all 256 selectors 0.
If they are right, my `A` is meaningless and the four beating pairs are roots of the wrong system.

Three checks against the ONE verified object in this lab, the 39,026 deliverable:
  (1) are any of the 253 pin/selector variables actually 1 in it?
  (2) do any ladder point coordinates appear as wire VALUES in it?
  (3) do any fold intermediates (partial sums of the chain) appear as wire values?
A model that cannot reproduce a verified point cannot be trusted to produce a new one.
"""
import sys, json, collections
sys.path.insert(0, '.')
import model
from model import P, S
lad = [(int(a), int(b)) for _, a, b in json.load(open('ladder.json'))['ladder']]
T = model.to_short(model.TARGET)
D = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
dv = {}
for k, v in D.items():
    i = int(k[2:]) if k.startswith('x_') else int(k)
    dv[i] = int(v)
pins = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/pins.json'))
PV = sorted(int(x) for x in pins)

print('=== (1) selector/pin variables set in the deliverable ===')
on = [v for v in PV if dv.get(v, 0) not in (0,)]
print('   pin variables total %d ; nonzero in deliverable: %d  %s' % (len(PV), len(on), on[:20]))
print('   x24601 = %s   x2081 = %s   (my claimed ON-set)' % (dv.get(24601, 0), dv.get(2081, 0)))

print('\n=== (2) do ladder coordinates appear as wire values? ===')
vals = collections.Counter(dv.values())
valset = set(dv.values())
# wires may hold the unshifted coordinate (short - S) or the short one
hits = []
for i, (x, y) in enumerate(lad):
    for cand, tag in (((x - S) % P, 'x-unshifted'), (x, 'x-short'), (y, 'y')):
        if cand in valset: hits.append((i, tag))
print('   ladder points with a coordinate present on some wire: %d / 256' % len({h[0] for h in hits}))
print('   sample:', hits[:10])
print('   target coords present?  x:%s  y:%s'
      % (((T[0] - S) % P) in valset or T[0] in valset, T[1] in valset))

print('\n=== (3) do chain fold intermediates appear as wire values? ===')
def chord(A, B):
    ax, ay = A; bx, by = B
    if (ax - bx) % P == 0: return None
    l = (by - ay) * pow(bx - ax, P - 2, P) % P
    sx = (l * l - ax - bx) % P
    return (sx, (l * (ax - sx) - ay) % P)
acc = lad[0]; found = 0; tested = 0
for i in (72, 235):
    nxt = chord(acc, lad[i])
    if nxt is None: continue
    tested += 1
    if nxt[0] in valset or ((nxt[0] - S) % P) in valset or nxt[1] in valset: found += 1
    acc = nxt
print('   fold intermediates for my ON-set {72,235}: %d of %d appear on a wire' % (found, tested))
print('   final fold value on a wire? %s' % (acc[0] in valset or ((acc[0] - S) % P) in valset))
json.dump({'pins_nonzero_in_deliverable': on, 'x24601': dv.get(24601, 0), 'x2081': dv.get(2081, 0),
           'ladder_coord_wire_hits': len({h[0] for h in hits}),
           'fold_intermediates_found': found}, open('runs/validate_A.json', 'w'), indent=1)
