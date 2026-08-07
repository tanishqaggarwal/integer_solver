"""S10 step 87: the ANCESTOR CONE of the two binding residues.

D0 = (x_7068 - x_2099) mod p and K2 = x_28730 mod p are the whole problem.
Walk back through `definer` to every variable that can influence them, and
inside that cone find every LOAD PIN -- an atom that ties a variable to a
constant -- and in particular the boolean-GATED pins b*(x - CONST), whose b=0
branch would unpin a constant and hand us a unit-granularity free input.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
P = ad.P
definer, atom_out = L.definer, L.atom_out

BOOL = set()
for _a, _poly in enumerate(L.polys):
    _ks = list(_poly.items())
    if len(_ks) == 2:
        _sq = [m for m, c in _ks if len(m) == 2 and m[0] == m[1]]
        _li = [m for m, c in _ks if len(m) == 1]
        if _sq and _li and _sq[0][0] == _li[0][0]:
            BOOL.add(_li[0][0])

TARGETS = [7068, 2099, 28730]
for t in TARGETS:
    print(f'x_{t}: {"GATE OUTPUT of a"+str(definer[t]) if t in definer else "FREE INPUT"}')

# ---- ancestor cone -------------------------------------------------------
cone = set()
stack = list(TARGETS)
while stack:
    t = stack.pop()
    if t in cone: continue
    cone.add(t)
    a = definer.get(t)
    if a is None: continue
    for w in L.avars[a]:
        if w != t and w not in cone:
            stack.append(w)
print(f'\nancestor cone of the three targets: {len(cone)} variables '
      f'({sum(1 for t in cone if t in definer)} gate outputs, '
      f'{sum(1 for t in cone if t not in definer)} free inputs)')
print(f'   booleans inside the cone: {len(cone & BOOL)}')

# ---- every atom fully inside the cone ------------------------------------
inside = [a for a in range(L.NA) if L.avars[a] and set(L.avars[a]) <= cone]
print(f'atoms wholly inside the cone: {len(inside)}')

# ---- classify pins -------------------------------------------------------
def big_consts(a):
    return sorted({abs(c) for c in L.polys[a].values() if abs(c) > 2**40})

gated, bare, other = [], [], []
for a in inside:
    poly = L.polys[a]
    vs = set(L.avars[a])
    bs = vs & BOOL
    if len(poly) <= 3 and bs:
        # b*(x - CONST) shape: monomials (b,x) and (b,)
        ms = list(poly)
        if all(len(m) <= 2 for m in ms) and any(len(m) == 2 for m in ms):
            gated.append(a); continue
    if len(vs) == 1 and len(poly) <= 2 and big_consts(a):
        bare.append(a); continue
    other.append(a)
print(f'   boolean-GATED pins inside the cone : {len(gated)}')
print(f'   bare constant pins inside the cone : {len(bare)}')

def show(a):
    return T.fmt(a) if hasattr(T, 'fmt') else str(L.polys[a])[:110]

print('\n--- boolean-gated pins in the cone (first 25) ---')
for a in gated[:25]:
    bs = sorted(set(L.avars[a]) & BOOL)
    print(f'  a{a}  vars={sorted(set(L.avars[a]))}  bool={bs}  eqs={len(L.atom2eq[a])}')
print('\n--- bare constant pins in the cone (first 25) ---')
for a in bare[:25]:
    print(f'  a{a}  var={sorted(set(L.avars[a]))}  eqs={len(L.atom2eq[a])}  '
          f'const~{[str(c)[:20] for c in big_consts(a)]}')
json.dump({'cone': sorted(cone), 'gated': gated, 'bare': bare},
          open(os.path.join(HERE, 'cone.json'), 'w'))
print('\nsaved cone.json')
