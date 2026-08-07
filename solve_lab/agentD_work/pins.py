"""Census of atom shapes: constant pins, boolean checks, wire pins, handles."""
import collections, sys
import dlib as L

P = L.P

shapes = collections.Counter()
pins = []       # (atom, var, const)
for a in range(L.NA):
    p = L.polys[a]
    deg = max((len(m) for m in p), default=0)
    nt = len(p)
    lin = [m for m in p if len(m) == 1]
    const = p.get((), 0)
    if deg == 1 and len(lin) == 1 and const:
        pins.append((a, lin[0][0], p[lin[0]], const))
    shapes[(deg, nt)] += 1

print('total atoms', L.NA)
print('constant pins (c*x + K):', len(pins))
big = [(a, u, c, k) for (a, u, c, k) in pins if abs(k) > 10**6]
print('  with |K| > 1e6:', len(big))
for a, u, c, k in sorted(big, key=lambda t: -abs(t[3]))[:40]:
    print(f'   a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} x_{u} = {-k//c if k % c == 0 else "?"}   bits={abs(k).bit_length()}')
print()
sm = [(a, u, c, k) for (a, u, c, k) in pins if abs(k) <= 10**6]
cnt = collections.Counter((-k // c if k % c == 0 else None) for a, u, c, k in sm)
print('small-pin value histogram:', cnt.most_common(12))
