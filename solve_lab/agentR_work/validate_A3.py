#!/usr/bin/env python3
"""Checks (2)/(3) redone: wire values are UNREDUCED (~89 digits), so searching for the
reduced/shifted coordinate was a false negative -- my own bug, same class as the pins lookup.
Reduce every wire value mod P and search in both frames.

Real question for solve2.py: is there an ACCUMULATOR CHAIN, and what is it seeded at?
The deliverable has exactly leaves 72 and 235 live, and leaf 0 is NOT live -- so a chain seeded
at L_0 (what solve2.py assumes) would be wrong.  Test what actually appears on wires."""
import sys, json
sys.path.insert(0, '.')
import model
from model import P, S
D = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
dv = {}
for k, v in D.items():
    dv[int(k[2:]) if k.startswith('x_') else int(k)] = int(v)
red = {}
for w, v in dv.items(): red.setdefault(v % P, []).append(w)
lad = [(int(a), int(b)) for _, a, b in json.load(open('ladder.json'))['ladder']]
T = model.to_short(model.TARGET)
def hit(val):
    for c in (val % P, (val - S) % P, (val + S) % P):
        if c in red: return red[c][:3]
    return None
def chord(A, B):
    ax, ay = A; bx, by = B
    if (ax - bx) % P == 0: return None
    l = (by - ay) * pow(bx - ax, P - 2, P) % P
    sx = (l * l - ax - bx) % P
    return (sx, (l * (ax - sx) - ay) % P)

print('distinct reduced wire values: %d' % len(red))
n = sum(1 for i, (x, y) in enumerate(lad) if hit(x) or hit(y))
print('ladder points with a coordinate on a wire (reduced): %d / 256' % n)
for i in (0, 72, 235):
    print('   leaf %-4d x-> %s   y-> %s' % (i, hit(lad[i][0]), hit(lad[i][1])))
print('target  x-> %s   y-> %s' % (hit(T[0]), hit(T[1])))

print('\nTHE TEST: is the fold of the two live leaves on a wire?')
Ssum = chord(lad[72], lad[235])
print('   L72 + L235          x-> %s   y-> %s' % (hit(Ssum[0]), hit(Ssum[1])))
c0 = chord(lad[0], lad[72])
print('   L0 + L72 (my seed)  x-> %s   y-> %s' % (hit(c0[0]), hit(c0[1])))
if c0:
    c1 = chord(c0, lad[235])
    print('   L0+L72+L235         x-> %s   y-> %s' % (hit(c1[0]), hit(c1[1])))
json.dump({'ladder_on_wire': n, 'fold_two_live': bool(hit(Ssum[0])),
           'fold_with_L0_seed': bool(hit(c1[0])) if c0 else None}, open('runs/validate_A3.json', 'w'), indent=1)
