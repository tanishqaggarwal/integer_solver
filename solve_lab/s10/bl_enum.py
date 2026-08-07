"""bl_enum: enumerate boolean free inputs, MUX controls, and load pins."""
import os, sys, json, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

# ---------- 1. booleans ----------
BOOL = set()          # var -> its boolean atom
BOOLATOM = {}
for a, p in enumerate(L.polys):
    ks = list(p.items())
    if len(ks) == 2:
        sq = [m for m, c in ks if len(m) == 2 and m[0] == m[1]]
        li = [m for m, c in ks if len(m) == 1]
        if sq and li and sq[0][0] == li[0][0]:
            BOOL.add(li[0][0]); BOOLATOM[li[0][0]] = a
print(f'boolean variables: {len(BOOL)}   free ones: {len(BOOL & FREE)}')

# ---------- 2. ancestor cones ----------
RES = [21617, 29539]
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]

def cone_of_atoms(atoms):
    c, st = set(), []
    for a in atoms: st += list(L.avars[a])
    while st:
        t = st.pop()
        if t in c: continue
        c.add(t)
        d = definer.get(t)
        if d is None: continue
        for w in L.avars[d]:
            if w != t: st.append(w)
    return c

C_res = cone_of_atoms(RES)
C_seven = cone_of_atoms(SEVEN)
C_all = C_res | C_seven
print(f'cone(a21617,a29539) vars {len(C_res)};  cone(seven) {len(C_seven)};  union {len(C_all)}')

B_res = sorted(C_res & BOOL & FREE)
B_seven = sorted(C_seven & BOOL & FREE)
B_all = sorted(C_all & BOOL & FREE)
print(f'boolean free inputs: in res-cone {len(B_res)}, in seven-cone {len(B_seven)}, union {len(B_all)}')

# ---------- 3. MUX-ness: does the boolean appear in a PRODUCT monomial? ----------
prod_use = collections.defaultdict(list)   # var -> [(atom, partner_vars)]
for a, p in enumerate(L.polys):
    if a in BOOLATOM.values() and len(p) == 2: pass
    for m, c in p.items():
        if len(m) >= 2:
            s = set(m)
            for u in s:
                if u in BOOL:
                    prod_use[u].append((a, tuple(sorted(x for x in m if x != u))))

# ---------- 4. LOAD PINS: atom == b*x - K*b  (2 monomials: (b,x) and (b,)) ----------
# more generally: every monomial contains b, and after factoring b the poly is
# linear in exactly one other var with a constant term.
PINS = []   # (atom, b, x, const)
for a, p in enumerate(L.polys):
    ms = list(p.items())
    if len(ms) != 2: continue
    # find pattern c1*(b,x) + c0*(b,)
    m2 = [(m, c) for m, c in ms if len(m) == 2]
    m1 = [(m, c) for m, c in ms if len(m) == 1]
    if len(m2) != 1 or len(m1) != 1: continue
    (mm, c1), = m2
    (nn, c0), = m1
    b = nn[0]
    if b not in mm: continue
    if mm[0] == mm[1]: continue      # that's the boolean atom itself
    x = mm[0] if mm[1] == b else mm[1]
    if x == b: continue
    if c0 % c1: continue
    PINS.append((a, b, x, -c0 // c1))
print(f'\nload pins of the form b*(x - K): {len(PINS)}')
pins_bool = [t for t in PINS if t[1] in BOOL]
print(f'   with a BOOLEAN gate variable b: {len(pins_bool)}')
pins_bool_free = [t for t in pins_bool if t[1] in FREE]
print(f'   with b a FREE boolean:          {len(pins_bool_free)}')
pins_in_cone = [t for t in pins_bool if t[3] != 0 and (t[2] in C_all or t[1] in C_all)]
print(f'   b or x inside the union cone:   {len(pins_in_cone)}')

json.dump({
  'BOOL': sorted(BOOL), 'BOOLFREE': sorted(BOOL & FREE),
  'B_res': B_res, 'B_seven': B_seven, 'B_all': B_all,
  'C_res': sorted(C_res), 'C_seven': sorted(C_seven),
  'PINS': [[a, b, x, str(k)] for a, b, x, k in PINS],
}, open(os.path.join(HERE, 'bl_enum.json'), 'w'))

# ---------- report ----------
print('\n--- boolean free inputs in the union cone, MUX-ness ---')
rows = []
for u in B_all:
    pu = prod_use.get(u, [])
    # exclude the boolean's own atom b*b - b
    pu = [(a, prt) for a, prt in pu if a != BOOLATOM.get(u)]
    npin = sum(1 for (aa, bb, xx, kk) in PINS if bb == u)
    rows.append((len(pu), npin, u, pu[:3]))
rows.sort(reverse=True)
print(f'{"var":>8} {"#prodUses":>9} {"#pins":>6}  sample')
for n, npin, u, s in rows[:40]:
    print(f'x_{u:<7} {n:>9} {npin:>6}  {[(a, prt) for a, prt in s]}')
nmux = sum(1 for n, npin, u, s in rows if n)
print(f'\nof {len(B_all)} cone booleans, {nmux} appear in a product (MUX-like), '
      f'{len(B_all)-nmux} appear only linearly')
