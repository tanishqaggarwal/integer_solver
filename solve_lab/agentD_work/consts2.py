"""Full census of large literals in atoms (constant terms AND big coefficients)."""
import collections
import dlib as L
P = L.P

rows = []
for a in range(L.NA):
    p = L.polys[a]
    big = {m: c for m, c in p.items() if abs(c) > 10**20}
    if big:
        rows.append((a, big))
print('atoms containing a literal > 1e20:', len(rows))
vals = collections.Counter()
for a, big in rows:
    for m, c in big.items():
        vals[abs(c)] += 1
print('distinct big literals:', len(vals))
for k, n in sorted(vals.items(), key=lambda t: -t[1])[:40]:
    print(f'  n={n:<4} bits={k.bit_length():<4} k%p={k % P}  {k}')
print()
print('--- per-atom (first 60, sorted by #equations desc) ---')
rows.sort(key=lambda t: -len(L.atom2eq.get(t[0], {})))
for a, big in rows[:60]:
    print(f'a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} gate={a in L.atom_out} :: {L.atom_src[a][:200]}')
