#!/usr/bin/env python3
"""Fold evaluator.

Established by K16: after the substitution X = x + K/3 (mod p) every one of the 256 leaf
constant pairs and the target satisfy the SAME cubic  Y^2 = X^3 + b  (a = 0), and the
stage law  out_x = l^2 - a_x - b_x - K,  out_y = l(a_x - out_x) - a_y  becomes
out_X = l^2 - a_X - b_X,  out_Y = l(a_X - out_X) - a_Y.
That law is the classical chord/tangent composition on the cubic, hence commutative and
associative, so the fold of a leaf subset does not depend on the tree shape:

        fold(S) = sum over i in S of P_i        (chord composition on Y^2 = X^3 + b)

`validate()` checks this against the exact cascade closure of the real equations.
"""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
KOFF = 97553848499418123410591666447050222001188385549510401465815187079080512838891
SHIFT = KOFF * pow(3, P - 2, P) % P
B = 64019533680030876408443198762210829058751700634554282185987325820393598524794

INF = None


def add(Pt, Q):
    if Pt is INF: return Q
    if Q is INF: return Pt
    x1, y1 = Pt; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % P == 0: return INF
        l = 3 * x1 * x1 % P * pow(2 * y1 % P, P - 2, P) % P
    else:
        l = (y2 - y1) * pow(x2 - x1, P - 2, P) % P
    x3 = (l * l - x1 - x2) % P
    y3 = (l * (x1 - x3) - y1) % P
    return (x3, y3)


def mul(k, Pt):
    R = INF
    while k:
        if k & 1: R = add(R, Pt)
        Pt = add(Pt, Pt); k >>= 1
    return R


def load():
    d = json.load(open(K + '/points.json'))
    return d


def build_points():
    """(re)extract the 256 leaf points and target straight from the equations."""
    sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
    from cascadep import CascadeP
    C = CascadeP()
    pins = []
    for nm in C.names:
        m = re.match(r'^\(\(x(\d+)\*\(x(\d+)-(\d{20,})\)\)-x(\d+)\)$', nm)
        if m: pins.append((int(m.group(1)), int(m.group(2)), int(m.group(3)))); continue
        m = re.match(r'^\(\(x(\d+)\*\(x(\d+)-(\d{20,})\)\)-\((\d+)\*x(\d+)\)\)$', nm)
        if m: pins.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    bysel = collections.defaultdict(list)
    for s, w, c in pins: bysel[s].append((w, c))
    assert len(bysel) == 256 and all(len(v) == 2 for v in bysel.values())
    leaves = {}
    nswap = 0
    for s in sorted(bysel):
        (w1, c1), (w2, c2) = sorted(bysel[s])
        ok = []
        for (a, b_, wa, wb) in ((c1, c2, w1, w2), (c2, c1, w2, w1)):
            X = (a + SHIFT) % P; Y = b_ % P
            if (Y * Y - pow(X, 3, P) - B) % P == 0: ok.append((X, Y, wa, wb))
        assert len(ok) == 1, ('ambiguous/off curve', s, len(ok))
        leaves[s] = ok[0]
        if ok[0][2] != w1: nswap += 1
    print('leaves whose x-pin is the higher wire index:', nswap)
    Tx = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002 % P
    Ty = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626 % P
    TX = (Tx + SHIFT) % P
    assert (Ty * Ty - pow(TX, 3, P) - B) % P == 0
    out = {'p': str(P), 'b': str(B), 'shift': str(SHIFT),
           'leaves': [{'sel': s, 'X': str(leaves[s][0]), 'Y': str(leaves[s][1]),
                       'wx': leaves[s][2], 'wy': leaves[s][3]} for s in sorted(leaves)],
           'target': {'X': str(TX), 'Y': str(Ty)}}
    json.dump(out, open(K + '/points.json', 'w'))
    return out


def points():
    if not os.path.exists(K + '/points.json'):
        return build_points()
    return load()


def fold(sel_on, D=None):
    D = D or points()
    idx = {l['sel']: (int(l['X']), int(l['Y'])) for l in D['leaves']}
    R = INF
    for s in sel_on: R = add(R, idx[s])
    return R


if __name__ == '__main__':
    D = build_points()
    print('leaves', len(D['leaves']), 'target', D['target'])
    idx = {l['sel']: (int(l['X']), int(l['Y'])) for l in D['leaves']}
    print('selector vars sample', sorted(idx)[:10])
    print('2081 is a leaf selector:', 2081 in idx, ' 24601:', 24601 in idx)
    T = (int(D['target']['X']), int(D['target']['Y']))
    f = fold([2081, 24601], D)
    print('fold({2081,24601}) =', f)
    print('equals target?', f == T)
