#!/usr/bin/env python3
"""K22: carry DP for a root degeneracy.

Leaf e has point 2^e * G (chain from K19).  The root's A half owns exponent set IA,
the B half owns IB, and IA + IB partition {0..255}.  A degenerate root needs
    sum_{e in SA} 2^e  ==  sum_{f in SB} 2^f   (mod N)
with SA subset IA, SB subset IB, both non-empty.  Writing x = sum_SA 2^e and
y = sum_SB 2^f we have 0 <= x,y < 2^256 and N > 2^255, so the congruence forces
x - y = +-N exactly.  Each bit position belongs to exactly one side, so the schoolbook
addition x = y + N leaves NO freedom: the DP is deterministic and just has to be run.
"""
import sys, os, json
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
import fold as FD

N = int(json.load(open(K + '/order.json'))['N'][0])
rs = json.load(open(K + '/rootsupport.json'))
ch = json.load(open(K + '/chain.json'))
sel2exp = {}
for i in range(256): sel2exp[ch['sel'][str(i)]] = ch['exp'][str(i)]
IA = set(sel2exp[s] for s in rs['A.x']) | set(sel2exp[s] for s in rs['A.y']) | {163}
IB = set(sel2exp[s] for s in rs['B.x']) | set(sel2exp[s] for s in rs['B.y'])
# exponent 163 measured A-side in K26 by driving the real equations with ON={x28785,x2081}
print('|IA| =', len(IA), '|IB| =', len(IB), 'overlap', len(IA & IB), 'union', len(IA | IB))
assert not (IA & IB) and len(IA | IB) == 256


def solve(m, plusset, minusset):
    """find x with bits in plusset, y with bits in minusset, x = y + m.  Deterministic."""
    x = 0; y = 0; c = 0
    for i in range(256):
        mi = (m >> i) & 1
        if i in plusset:
            s = mi + c            # y_i = 0
            x |= (s & 1) << i
            c = s >> 1
        else:
            s = mi + c            # x_i = 0, choose y_i so that s+y_i is even
            if s == 0: yi, c = 0, 0
            elif s == 1: yi, c = 1, 1
            else: yi, c = 0, 1
            y |= yi << i
    if c: return None
    assert x - y == m
    return x, y


out = {}
for lab, m, plus, minus in [('x-y=+N', N, IA, IB), ('y-x=+N', N, IB, IA)]:
    r = solve(m, plus, minus)
    print(lab, '->', 'no solution (final carry)' if r is None else 'SOLVED')
    if r:
        x, y = r
        if lab == 'x-y=+N': X, Y = x, y
        else: X, Y = y, x
        print('   x (A-side multiplier) =', X)
        print('   y (B-side multiplier) =', Y)
        print('   x mod N ==', X % N, '  y mod N ==', Y % N, ' equal:', X % N == Y % N)
        print('   popcount A =', bin(X).count('1'), ' popcount B =', bin(Y).count('1'))
        out[lab] = {'A_mult': str(X), 'B_mult': str(Y),
                    'SA_exp': [i for i in range(256) if (X >> i) & 1],
                    'SB_exp': [i for i in range(256) if (Y >> i) & 1]}
        # verify with the actual points
        D = FD.points()
        bypow = {}
        for i_s, e in ch['exp'].items():
            bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
        fa = FD.INF
        for e in out[lab]['SA_exp']: fa = FD.add(fa, bypow[e])
        fb = FD.INF
        for e in out[lab]['SB_exp']: fb = FD.add(fb, bypow[e])
        print('   fold(SA) =', fa)
        print('   fold(SB) =', fb)
        print('   *** EQUAL:', fa == fb, ' non-identity:', fa is not None)
json.dump(out, open(K + '/degenerate.json', 'w'))
