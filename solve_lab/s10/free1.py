"""S11 step 43: can congruence 1 be freed?

A0 + 7376877*A6 == C0 (mod p) exists only because x_7068 is pinned mod p, and the
ONLY thing moving x_7068 breaks is a29539 (13 equations, measured).  If a29539's
congruence can be re-fixed for less than 1 equation, c drops to 1 and
failing = 12 - 7 + 1 = 6  ->  39,027.  Price every repair in a29539's support.
"""
import os, sys, time
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
print(f'frame 2 base {score(base)}')
v = list(base); v[7068] = v[7068] + 1          # move C0 off its residue
fwd(v, rounds=8)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a] and a not in SSET]
print(f'x_7068 += 1 -> broken outside the seven: {nz}')
vm = [x % P for x in v]
g = grad(29539, vm)
print(f"a29539's gradient support: {len(g)} free inputs", flush=True)
r = av[29539] % P
def close_handle(w):
    """absorb a29539/p through x_30163 once its congruence holds."""
    t = L.atom_out[29538][1]
    tgt = T.solve_lin(29539, t, w)
    if tgt is None: return None
    vv = list(w); vv[t] = tgt
    nv = T.solve_lin(29538, 30163, vv)
    if nv is None: return None
    z = list(w); z[30163] = nv
    fwd(z, rounds=8)
    return z
res = []
t0 = time.time()
for i, (u, d) in enumerate(sorted(g.items(), key=lambda kv: len(L.var_atoms[kv[0]]))):
    if u in (2081, 4287) or d % P == 0: continue
    delta = (-r * pow(d, -1, P)) % P
    w = list(v); w[u] = w[u] + delta
    fwd(w, rounds=8)
    aw = L.all_atom_values(w)
    if aw[29539] % P: continue
    z = close_handle(w)
    if z is None: z = w
    az = L.all_atom_values(z)
    out = [a for a in range(L.NA) if az[a] and a not in SSET]
    eqs = set()
    for a in out: eqs |= set(L.atom2eq[a])
    res.append((len(eqs), u, out))
    if len(eqs) <= 13:
        print(f'  x_{u:<7} (cons {len(L.var_atoms[u]):>2}) -> outside-seven {out} '
              f'({len(eqs)} eqs)  score {score(z)}', flush=True)
    if len(eqs) == 0:
        print('    *** a29539 RE-FIXED AT ZERO COST -- congruence 1 is free')
        T.save(z, os.path.join(HERE, 'free1.json'))
    if i > 90: break
res.sort()
print(f'\n{len(res)} repairs priced ({time.time()-t0:.0f}s); cheapest:')
for k, u, out in res[:8]:
    print(f'  {k:>4} equations  via x_{u}  {out}')
