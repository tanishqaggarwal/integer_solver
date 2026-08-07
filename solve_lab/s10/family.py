"""S10 step 105: enumerate the WHOLE 'A == B (mod p)' gadget family and its graph."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
av = L.all_atom_values(v)

fam = []
for a in range(L.NA):
    poly = L.polys[a]
    if len(poly) != 3: continue
    if any(len(m) != 1 for m in poly): continue
    items = sorted(poly.items(), key=lambda kv: -abs(kv[1]))
    (m1, c1), (m2, c2), (m3, c3) = items
    if c1 == -c2 and abs(c3) == 1 and abs(c1) > 1:
        fam.append((a, m1[0], m2[0], m3[0], c1))
print(f'linear gadgets  c*(A-B) - C : {len(fam)}')
print(f'   currently failing: {[a for a,_,_,_,_ in fam if av[a]]}')
# is C always p*handle?
pq = 0
for a, A, B, C, c in fam:
    d = definer.get(C)
    if d is None: continue
    vs = [w for w in set(L.avars[d]) if w != C]
    if len(vs) == 2 and any(v[w] == P for w in vs):
        h = [w for w in vs if v[w] != P]
        if h and len(L.var_atoms[h[0]]) == 1: pq += 1
print(f'   with C = p * (solo free handle): {pq} / {len(fam)}')

# graph: gadget -> gadget, via shared ancestor free inputs
vm = [x % P for x in v]
supp = {}
for a, A, B, C, c in fam:
    supp[a] = set(ad.grad(a, vm))
print(f'\ngadget -> which other gadgets share a moving free input')
ids = [a for a, *_ in fam]
deg = collections.Counter()
for a in ids:
    nb = [b for b in ids if b != a and supp[a] & supp[b]]
    deg[a] = len(nb)
print(f'   degree distribution: {collections.Counter(deg.values())}')
for a, A, B, C, c in fam:
    nb = [b for b in ids if b != a and supp[a] & supp[b]]
    print(f'  a{a:<6} ({len(L.atom2eq[a]):>2} eqs) {c}*(x_{A} - x_{B}) - x_{C}   '
          f'{"FAIL" if av[a] else "ok":>4}  |grad|={len(supp[a]):<4} nbrs={nb}')
