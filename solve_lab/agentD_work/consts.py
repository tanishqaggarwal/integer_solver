"""All atoms containing a large integer constant, plus their shape."""
import collections
import dlib as L
P = L.P

rows = []
for a in range(L.NA):
    p = L.polys[a]
    big = [c for m, c in p.items() if abs(c) > 10**12]
    k = p.get((), 0)
    if abs(k) > 10**12 or big:
        rows.append((a, k, big))

print('atoms with a >1e12 coefficient or constant:', len(rows))
constrows = [(a, k) for a, k, b in rows if abs(k) > 10**12]
print('  with big CONSTANT term:', len(constrows))
for a, k in sorted(constrows, key=lambda t: -abs(t[1])):
    print(f'  a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} K={k}  bits={abs(k).bit_length()}  K%p={k % P}')
    print(f'        {L.atom_src[a][:220]}')
