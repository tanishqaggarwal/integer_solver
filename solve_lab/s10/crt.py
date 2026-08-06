"""S10 step 98: CRT-corrected Newton step.

The Newton shift delta is pinned only mod p, so x_6418 += delta + k*p is free for
every k.  Choose k so that the SECOND condition -- 15804267 | (x_6418 - C), which
a3576 needs in order to close through x_26777 -- also holds.  One congruence mod p,
one mod 15804267, gcd 1: CRT gives both at once.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from newton_modp import newton_moves, score, FORBID
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
print(f'start score {score(v0)}')
print(f'x_26777 free? {26777 in FREE}  consumers {len(L.var_atoms[26777])}: '
      f'{sorted(L.var_atoms[26777])}')

# a3576 = -C*x_2081 - 15804267*x_26777 + x_2081*x_6418
C = -[c for m, c in L.polys[3576].items() if m == (2081,)][0]
M = 15804267
print(f'C has {C.bit_length()} bits;  M = {M}')

av0 = L.all_atom_values(v0)
delta = dict(newton_moves(29539, v0, av0, topn=200))[6418]
base = v0[6418] + delta
k = (-(base - C) * pow(P, -1, M)) % M
print(f'delta mod p ok; k = {k} so that M | (x_6418 - C)')
v = list(v0); v[6418] = base + k * P
assert (v[6418] - C) % M == 0
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
print(f'after x_6418 move: score {score(v)}  nonzero {[a for a in range(L.NA) if av[a]]}')
print(f'   a29539 mod p == 0 ? {av[29539] % P == 0}')

# close a29539 through its p-handle x_30163
w = atom_out[29538][1]
tgt = T.solve_lin(29539, w, v)
nv = T.solve_lin(29538, 30163, [x if i != w else tgt for i, x in enumerate(v)])
print(f'   handle x_30163 -> {str(nv)[:40]}')
v[30163] = nv
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
print(f'after handle: score {score(v)}  nonzero {[a for a in range(L.NA) if av[a]]}')

# close a3576 through x_26777
nv = T.solve_lin(3576, 26777, v)
print(f'   x_26777 -> {str(nv)[:40]}  (solvable: {nv is not None})')
if nv is not None:
    v[26777] = nv
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    s = score(v)
    print(f'after x_26777: score {s}  nonzero {[a for a in range(L.NA) if av[a]]}')
    T.save(v, os.path.join(HERE, f'crt_{s}.json'))
    print(f'saved crt_{s}.json')
