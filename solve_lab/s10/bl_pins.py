"""bl_pins: general search for boolean-gated pins (atoms divisible by a boolean b)."""
import os, sys, json, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
E = json.load(open(os.path.join(HERE, 'bl_enum.json')))
BOOL = set(E['BOOL']); BOOLFREE = set(E['BOOLFREE'])
B_all = E['B_all']; C_res = set(E['C_res']); C_seven = set(E['C_seven'])
CONE = C_res | C_seven

# --- atoms divisible by a boolean variable b (every monomial contains b) ---
DIV = collections.defaultdict(list)     # b -> [atom]
for a, p in enumerate(L.polys):
    ms = list(p.keys())
    if len(ms) < 2: continue
    common = set(ms[0])
    for m in ms[1:]:
        common &= set(m)
        if not common: break
    for b in common:
        if b in BOOL: DIV[b].append(a)
print(f'booleans that divide some atom: {len(DIV)}')

# For each such atom, quotient by b and inspect: is it (x - K)?
def quotient(p, b):
    q = collections.defaultdict(int)
    for m, c in p.items():
        mm = list(m); mm.remove(b)
        q[tuple(sorted(mm))] += c
    return dict(q)

pinlike = []      # (atom, b, x, K)  quotient = c*(x - K)
gatelike = []     # atom, b, quotient with >1 term
for b, As in DIV.items():
    for a in As:
        if len(L.polys[a]) == 2 and a == E and False: continue
        q = quotient(L.polys[a], b)
        ks = list(q.items())
        if len(ks) == 2:
            lin = [(m, c) for m, c in ks if len(m) == 1]
            con = [(m, c) for m, c in ks if len(m) == 0]
            if len(lin) == 1 and len(con) == 1:
                (mm, c1), = lin; (_, c0), = con
                if c0 % c1 == 0:
                    pinlike.append((a, b, mm[0], -c0 // c1)); continue
        gatelike.append((a, b, len(ks)))
print(f'  pin-like  b*(c1)*(x - K): {len(pinlike)}')
print(f'  other b*Q:               {len(gatelike)}')

sig = [t for t in pinlike if t[3] != 0]
print(f'  pins with NONZERO constant K: {len(sig)}')
for a, b, x, K in sig[:60]:
    inc = ('CONE' if (x in CONE or b in CONE) else '    ')
    fb = 'FREEB' if b in FREE else 'derB '
    dx = 'defX' if x in definer else 'freeX'
    print(f'  {inc} a{a:<6} b=x_{b:<6}({fb}) x=x_{x:<6}({dx}) K={str(K)[:60]}{"..." if len(str(K))>60 else ""}')

# --- unconditional constant pins: atom = c*x - c*K, 2 monomials, one is a bare const ---
CONST = []
for a, p in enumerate(L.polys):
    ks = list(p.items())
    if len(ks) != 2: continue
    lin = [(m, c) for m, c in ks if len(m) == 1]
    con = [(m, c) for m, c in ks if len(m) == 0]
    if len(lin) == 1 and len(con) == 1:
        (mm, c1), = lin; (_, c0), = con
        if c0 % c1 == 0: CONST.append((a, mm[0], -c0 // c1))
print(f'\nunconditional constant pins (x == K): {len(CONST)}')
inc = [t for t in CONST if t[1] in CONE]
print(f'  inside the union cone: {len(inc)}')
for a, x, K in inc[:40]:
    print(f'   a{a:<6} x_{x:<6} = {str(K)[:70]}')

json.dump({'pinlike': [[a, b, x, str(k)] for a, b, x, k in pinlike],
           'gatelike': gatelike,
           'const': [[a, x, str(k)] for a, x, k in CONST]},
          open(os.path.join(HERE, 'bl_pins.json'), 'w'))
