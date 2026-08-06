"""S11 step 2: the DEAD PATHS.

35,208 of 38,748 variables are 0 at our witness, so every monomial u*w with w = 0
has zero derivative in u.  The 1655x707 closure -- "full column rank, kernel 0" --
is therefore blind to every free input that reaches the cluster only through a
currently-dead product.  Find them: structural ancestors MINUS numeric gradient
support.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
BAD = [21617, 29539]

def struct_cone(seeds):
    cone, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in cone: continue
        cone.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
    return cone

seeds = set()
for a in BAD: seeds |= set(L.avars[a])
cone = struct_cone(seeds)
sfree = {u for u in cone if u in FREE}
nfree = set()
for a in BAD: nfree |= set(ad.grad(a, vm))
print(f'structural free ancestors of the cluster : {len(sfree)}')
print(f'numeric gradient support (nonzero d/du)  : {len(nfree & sfree)}')
dead = sorted(sfree - nfree)
print(f'DEAD free inputs (reach it, derivative 0): {len(dead)}')
print(f'   of those, currently zero: {sum(1 for u in dead if v[u] == 0)}')
print(f'   sample: {dead[:30]}')

# why is each dead?  which atom kills it
print('\nwhy they are dead (first 12):')
for u in dead[:12]:
    ats = sorted(L.var_atoms[u])
    reasons = []
    for a in ats:
        for m, c in L.polys[a].items():
            if u in m:
                others = [w for w in m if w != u]
                if others and all(v[w] == 0 for w in others):
                    reasons.append(f'a{a}: {u}*{others} with {others} = 0')
                    break
    print(f'  x_{u:<7} (val {v[u]}, in {len(ats)} atoms) {reasons[:2]}')

# the whole-instance version: how much of the circuit is switched off?
zero_vars = sum(1 for u in range(L.NVARS) if v[u] == 0)
dead_mon = live_mon = 0
for a in range(L.NA):
    for m, c in L.polys[a].items():
        if len(m) >= 2:
            if any(v[w] == 0 for w in m): dead_mon += 1
            else: live_mon += 1
print(f'\nquadratic monomials: {live_mon} live, {dead_mon} dead '
      f'({100*dead_mon/(live_mon+dead_mon):.1f}% switched off)')
