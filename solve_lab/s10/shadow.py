"""S10 step 100: are the 1-equation 'hardening' checks INDEPENDENT constraints,
or just Z-multiples of the gadget they shadow?"""
import os, sys, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'nonzero atoms and values:')
for a in nz:
    print(f'  a{a:<6} ({len(L.atom2eq[a])} eqs) {av[a]}')
print()
for s, g in [(37662, 21617), (40826, 29539)]:
    f = Fraction(av[s], av[g])
    print(f'a{s} / a{g} = {f}  -> integer multiple? {f.denominator == 1}')

# global census of the "A == B (mod p)" gadget family
fam = []
for a in range(L.NA):
    poly = L.polys[a]
    if len(poly) != 3: continue
    lins = [(m, c) for m, c in poly.items() if len(m) == 1]
    quads = [(m, c) for m, c in poly.items() if len(m) == 2]
    if len(lins) == 2 and len(quads) == 1 and lins[0][1] == -lins[1][1]:
        fam.append((a, lins, quads[0]))
print(f'\ngadget family  c*(A-B) - X*Y : {len(fam)} atoms')
sat = [a for a, _, _ in fam if av[a] == 0]
print(f'   satisfied at this state: {len(sat)} / {len(fam)}')
print(f'   FAILING: {[a for a, _, _ in fam if av[a]]}')
wire = collections.Counter()
handles = 0
for a, lins, (qm, qc) in fam:
    y, z = qm
    if v[y] == P or v[z] == P:
        wire[a] = 1
        h = z if v[y] == P else y
        if len(L.var_atoms[h]) == 1: handles += 1
print(f'   with a p-wire factor in the quadratic term: {sum(wire.values())}')
print(f'   ... of which the other factor is a SOLO handle: {handles}')
