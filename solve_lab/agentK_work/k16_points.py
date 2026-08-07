#!/usr/bin/env python3
"""K16: extract the 256 leaf points + target; test whether one cubic relation covers them all
(which makes the chord-with-offset law associative, i.e. the fold is a commutative group sum)."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from parse import node_str

KK = 97553848499418123410591666447050222001188385549510401465815187079080512838891
C = CascadeP()

pins = []   # (sel_var, wire_var, const)
pat1 = re.compile(r'^\(x(\d+)\*\(x(\d+)-(\d{20,})\)\)$')
for i, nm in enumerate(C.names):
    m = re.match(r'^\(\(x(\d+)\*\(x(\d+)-(\d{20,})\)\)-x(\d+)\)$', nm)
    if m:
        pins.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), 1)); continue
    m = re.match(r'^\(\(x(\d+)\*\(x(\d+)-(\d{20,})\)\)-\((\d+)\*x(\d+)\)\)$', nm)
    if m:
        pins.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(5)), int(m.group(4)))); continue
print('leaf pins:', len(pins))
bysel = collections.defaultdict(list)
for s, w, c, h, k in pins: bysel[s].append((w, c))
print('distinct selector vars:', len(bysel))
print('pins per selector histogram:', collections.Counter(len(v) for v in bysel.values()))

Tx = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002 % P
Ty = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626 % P
print('target', Tx, Ty)

# shift X = x + K/3 makes out_x = l^2 - a_x - b_x - K into out_X = l^2 - a_X - b_X
sh = KK * pow(3, P - 2, P) % P
pts = {}
for s, wl in bysel.items():
    if len(wl) != 2: continue
    (w1, c1), (w2, c2) = sorted(wl)
    pts[s] = ((w1, c1), (w2, c2))
print('leaves with a full pair:', len(pts))

# which of the two is x and which is y? try both orderings against a common cubic
def fitcurve(P1, P2, P3):
    """solve y^2 = X^3 + aX + b from two points, check third."""
    (X1, Y1), (X2, Y2) = P1, P2
    # Y1^2 - X1^3 = a X1 + b ; Y2^2 - X2^3 = a X2 + b
    r1 = (Y1 * Y1 - X1 ** 3) % P
    r2 = (Y2 * Y2 - X2 ** 3) % P
    if (X1 - X2) % P == 0: return None
    a = (r1 - r2) * pow(X1 - X2, P - 2, P) % P
    b = (r1 - a * X1) % P
    X3, Y3 = P3
    if (Y3 * Y3 - X3 ** 3 - a * X3 - b) % P != 0: return None
    return a, b

sels = sorted(pts)
cand = []
for order in (0, 1):
    L = []
    for s in sels:
        (w1, c1), (w2, c2) = pts[s]
        x, y = (c1, c2) if order == 0 else (c2, c1)
        L.append(((x + sh) % P, y % P))
    r = fitcurve(L[0], L[1], L[2])
    print('order', order, '->', r)
    if r:
        a, b = r
        bad = [i for i, (X, Y) in enumerate(L) if (Y * Y - X ** 3 - a * X - b) % P]
        print('   points off the curve:', len(bad), 'of', len(L))
        tgt = ((Tx + sh) % P, Ty % P)
        print('   target on curve:', (tgt[1] ** 2 - tgt[0] ** 3 - a * tgt[0] - b) % P == 0)
        if not bad:
            json.dump({'a': str(a), 'b': str(b), 'shift': str(sh), 'p': str(P),
                       'order': order,
                       'leaves': [{'sel': s, 'wx': pts[s][0][0] if order == 0 else pts[s][1][0],
                                   'x': str(L[i][0]), 'y': str(L[i][1])} for i, s in enumerate(sels)],
                       'target': {'X': str(tgt[0]), 'Y': str(tgt[1])}},
                      open(K + '/points.json', 'w'))
            print('   wrote points.json')
