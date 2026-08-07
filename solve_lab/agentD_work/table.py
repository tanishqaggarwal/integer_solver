"""Extract the 256 gated point-table entries and test them against secp256k1."""
import collections, sys, json
import dlib as L
P = L.P
A_, B_ = 0, 7

pins = []
for a in range(L.NA):
    p = L.polys[a]
    if max((len(m) for m in p), default=0) != 2:
        continue
    quad = {m: c for m, c in p.items() if len(m) == 2}
    lin = {m[0]: c for m, c in p.items() if len(m) == 1}
    if p.get(()) or len(quad) != 1:
        continue
    (qm, qc), = quad.items()
    if qc != 1 or qm[0] == qm[1]:
        continue
    bigl = [(u, c) for u, c in lin.items() if abs(c) > 10**20]
    if len(bigl) != 1:
        continue
    b, negC = bigl[0]
    if b not in qm:
        continue
    x = qm[0] if qm[1] == b else qm[1]
    pins.append((a, b, x, -negC))

bybit = collections.defaultdict(list)
for a, b, x, C in pins:
    bybit[b].append((a, x, C))
print('bits:', len(bybit))

oncurve = 0
tbl = {}
for b, lst in bybit.items():
    lst.sort()
    cs = [C % P for _, _, C in lst]
    tbl[b] = (lst, cs)
    # test both orderings for on-curve
    for xi, yi in ((0, 1), (1, 0)):
        X, Y = cs[xi], cs[yi]
        if (Y * Y - X * X * X - 7) % P == 0:
            oncurve += 1
            break
print('bit entries whose (C1,C2) is a secp256k1 point (either order):', oncurve, '/', len(bybit))

# report a few
for b, (lst, cs) in list(tbl.items())[:5]:
    X, Y = cs
    print(f'  bit x_{b}: X={X}  Y={Y}  onc={(Y*Y-X**3-7)%P==0} onc_sw={(X*X-Y**3-7)%P==0}')

json.dump({str(b): [[a, x, C] for a, x, C in lst] for b, (lst, cs) in tbl.items()},
          open('table.json', 'w'))
print('written table.json')
