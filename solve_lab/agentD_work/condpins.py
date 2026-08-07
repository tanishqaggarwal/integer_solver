"""Census of conditional pins  b*(x - C) - K*h  ."""
import collections, sys
import dlib as L
P = L.P

pins = []   # (atom, b, x, C, hvar/coeff)
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
    # need lin terms: -C*b  and  -K*h
    bigl = [(u, c) for u, c in lin.items() if abs(c) > 10**20]
    if len(bigl) != 1:
        continue
    b, negC = bigl[0]
    if b not in qm:
        continue
    x = qm[0] if qm[1] == b else qm[1]
    rest = [(u, c) for u, c in lin.items() if u != b]
    pins.append((a, b, x, -negC, rest))

print('conditional pins b*(x - C):', len(pins))
bs = collections.Counter(b for _, b, _, _, _ in pins)
print('distinct gating bits:', len(bs), ' pins/bit hist:', collections.Counter(bs.values()))
xs = collections.Counter(x for _, _, x, _, _ in pins)
print('distinct pinned vars:', len(xs), ' pins/var hist:', collections.Counter(xs.values()))
Cs = collections.Counter(C for _, _, _, C, _ in pins)
print('distinct constants:', len(Cs), ' bits/const hist:', collections.Counter(Cs.values()))
print()
freeb = [b for b in bs if b in L.freeset]
print('gating bits that are FREE inputs:', len(freeb))
print()
# where do pinned vars flow?
print('sample pins:')
for a, b, x, C, rest in pins[:10]:
    print(f'  a{a:<6} b=x_{b:<6}(free={b in L.freeset}) pins x_{x:<6} to {str(C)[:30]}... bits={C.bit_length()} rest={rest}')
