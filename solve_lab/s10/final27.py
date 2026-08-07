"""S11 step 20: A1 = 0 needs a7930's congruence fixed by something OTHER than
x_28730.  a7930's gradient support is tiny -- price every member."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad
P = ad.P
SSET = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
v0 = list(base); v0[28730] = v0[9413] * P        # A1 = 0
fwd(v0, rounds=8)
av0 = L.all_atom_values(v0)
print(f'A1 = 0; a22230 = {av0[22230]}')
vm = [x % P for x in v0]
g = grad(7930, vm)
print(f'a7930 gradient support: {sorted(g)}  (residue mod p '
      f'{"zero" if av0[7930] % P == 0 else "nonzero"})')

def handle_close(v):
    """close a7930 through x_7927's handle x_11052 if the congruence allows."""
    tgt = T.solve_lin(7930, 7927, v)
    if tgt is None: return None
    d = definer.get(7927)
    if d is None: return None
    vv = list(v); vv[7927] = tgt
    for u in sorted(set(L.avars[d])):
        if u == 7927 or u not in FREE: continue
        nv = T.solve_lin(d, u, vv)
        if nv is not None:
            w = list(v); w[u] = nv
            fwd(w, rounds=8)
            return w
    return None

r = av0[7930] % P
best = None
for u in sorted(g):
    if u == 28730: continue
    d = g[u] % P
    if d == 0: continue
    delta = (-r * pow(d, -1, P)) % P
    v = list(v0); v[u] = v[u] + delta
    fwd(v, rounds=8)
    av = L.all_atom_values(v)
    ok = av[7930] % P == 0
    w = handle_close(v) if ok else None
    if w is None:
        nz = [a for a in range(L.NA) if av[a] and a not in SSET]
    else:
        aw = L.all_atom_values(w)
        nz = [a for a in range(L.NA) if aw[a] and a not in SSET]
        v = w
    eqs = set()
    for a in nz: eqs |= set(L.atom2eq[a])
    av = L.all_atom_values(v)
    print(f'  x_{u:<7} (consumers {len(L.var_atoms[u])}) congruence fixed {ok}; '
          f'handle close {"OK" if w is not None else "no"};  '
          f'outside-seven {nz} ({len(eqs)} eqs);  a22230 '
          f'{"ZERO" if av[22230] == 0 else "nonzero"}', flush=True)
    if av[22230] == 0 and not nz:
        print('    *** A1 = 0 AND nothing else broken -- this is 39,027')
        T.save(v, os.path.join(HERE, 'A1zero.json'))
    if best is None or len(eqs) < best[0]:
        best = (len(eqs), u, v)
print(f'\ncheapest: x_{best[1]} leaving {best[0]} equations broken outside the seven')
