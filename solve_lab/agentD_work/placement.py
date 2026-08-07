"""Placement census: for every atom (and its equation set closed once), compute
  E   = equations touched
  n   = atoms entirely confined to E ("free" values)
  bound = |E| - n + c   (c = number of independent realisability congruences, >=1)
This is the score deficit lower bound for a residual placed there."""
import collections, sys, itertools
import dlib as L

atom_eqs = {a: set(L.atom2eq.get(a, {})) for a in range(L.NA)}
res = []
for a in range(L.NA):
    E = atom_eqs[a]
    if not E or len(E) > 30:
        continue
    atoms = set()
    for i in E:
        atoms |= set(L.eq_atoms[i][2])
    conf = [b for b in atoms if atom_eqs[b] <= E]
    res.append((len(E) - len(conf), len(E), len(conf), a))
res.sort()
print('atom            |E|  confined   |E|-n')
for d, ne, nc, a in res[:50]:
    print(f'  a{a:<6} {ne:<4} {nc:<4} {d:<4}  eqs={sorted(atom_eqs[a])[:6]}  :: {L.atom_src[a][:70]}')

print()
print('--- closure: grow E by unioning equations of confined atoms until stable ---')
def closure(a):
    E = set(atom_eqs[a])
    for _ in range(6):
        atoms = set()
        for i in E:
            atoms |= set(L.eq_atoms[i][2])
        conf = [b for b in atoms if atom_eqs[b] <= E]
        if not conf:
            break
        return E, conf
    return E, []

best = []
for d, ne, nc, a in res[:400]:
    E, conf = closure(a)
    best.append((len(E) - len(conf), len(E), len(conf), a))
best.sort()
for d, ne, nc, a in best[:30]:
    print(f'  a{a:<6} |E|={ne:<4} n={nc:<4} deficit>={d}')
