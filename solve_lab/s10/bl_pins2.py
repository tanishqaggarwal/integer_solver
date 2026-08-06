"""bl_pins2: direct enumeration of gated load pins  b*(x - K)  (b != x)."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, BOOLATOM, CANON, F2, pot, FORBID
P = 2**256-2**32-977

w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
v0 = L.load(os.path.join(HERE,'mod9118_0.json')); CANON.fwd(v0)

PINS = []
for a, p in enumerate(L.polys):
    ms = list(p.items())
    if len(ms) != 2: continue
    m2 = [(m, c) for m, c in ms if len(m) == 2]
    m1 = [(m, c) for m, c in ms if len(m) == 1]
    if len(m2) != 1 or len(m1) != 1: continue
    (mm, c1), = m2
    (nn, c0), = m1
    b = nn[0]
    if b not in mm: continue
    x = mm[0] if mm[1] == b else mm[1]
    if x == b: continue                       # boolean atom b*b-b
    if c0 % c1: continue
    PINS.append((a, b, x, -c0 // c1))
print(f'gated load pins b*(x - K), b != x : {len(PINS)}')
bb = [t for t in PINS if t[1] in BOOL]
print(f'   gate b is a BOOLEAN:        {len(bb)}')
bf = [t for t in bb if t[1] in CANON.FREE]
print(f'   gate b is a FREE boolean:   {len(bf)}')
print(f'   distinct gates: {sorted(set(t[1] for t in PINS))}')
print(f'   distinct pinned vars: {len(set(t[2] for t in PINS))}')
print()
for a, b, x, K in PINS:
    st = 'ACTIVE(pinning)' if w[b] == 1 else 'inactive'
    print(f'  a{a:<6} gate x_{b:<6}(={w[b]})  pins x_{x:<6} = {str(K)[:56]}...  bits={K.bit_length()}  [{st}]')
    print(f'        x_{x} currently = {str(w[x])[:56]}  bits={w[x].bit_length()}   K mod p = {K % P}')
    print(f'        x_{x} free? canon={x in CANON.FREE} f2={x in F2.FREE};  '
          f'#atoms using x_{x} = {len(L.var_atoms[x])}')

json.dump([[a, b, x, str(K)] for a, b, x, K in PINS], open(os.path.join(HERE,'bl_pins2.json'),'w'))

# ------ which pinned constants sit in the seven-check cone? ------
c2 = F2.cone([22229,22230,35758,35759,35760,35761,35762])
cc = CANON.cone([21617, 29539])
print('\npins whose pinned variable lies in a residual cone:')
for a, b, x, K in PINS:
    tags = []
    if x in c2: tags.append('SEVEN-CONE')
    if x in cc: tags.append('CLUSTER-CONE')
    if tags: print(f'  a{a} gate x_{b} pins x_{x}: {tags}')
