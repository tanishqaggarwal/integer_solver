"""S11 step 13: repair the cluster while holding the six-zero / a35758-nonzero
configuration.  If a7930 and a29539 both close, the score is 39,027."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import DETACH, definer, ORDER, FREE, CHECKS, fwd, grad
P = ad.P
FORBID = {2081, 4287}
PINNED = {29854, 31864, 10903, 1329, 642, 28730, 7068, 17325, 9413}

def renorm(v):
    v[8731] = (v[8731] // P) * P                 # keep p | x_8731
    v[31864] = -v[7075] * v[8731]
    v[10903] = v[31864] // P
    v[29854] = 5113045 * v[7075] * v[9118]
    v[1329] = v[29854] // P
    v[642] = v[17325] * P
    v[28730] = v[9413] * P
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    return v

def pot(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    return (L.NEQ - len(L.failing_eqs(av)), -len(nz)), av, nz

def moves(a, v, av, nnewton=24):
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID or w in PINNED: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in FREE: out.append((w, tgt))
        else:
            d = definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE or u in FORBID or u in PINNED: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: out.append((u, nv))
    r = av[a] % P
    if r:
        vm = [x % P for x in v]
        try: g = grad(a, vm)
        except Exception: g = {}
        c = sorted((len(L.var_atoms[u]), u, d) for u, d in g.items()
                   if u not in FORBID and u not in PINNED and d % P)
        for _, u, d in c[:nnewton]:
            step = (-r * pow(d, -1, P)) % P
            if u == 8731: step = (step // P) * P      # keep the congruence
            out.append((u, v[u] + step))
    return out

v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(v, rounds=6)
renorm(v)
cur, av, nz = pot(v)
print(f'start {cur[0]}  nonzero {nz}', flush=True)
t0 = time.time()
for it in range(40):
    got = None
    for a in nz:
        for u, nv in moves(a, v, av):
            tr = list(v); tr[u] = nv
            renorm(tr)
            p2, av2, nz2 = pot(tr)
            if p2 > cur: got = (a, u, p2, tr, av2, nz2); break
        if got: break
    if not got:
        print(f'it{it}: stuck at {cur[0]}  nonzero {nz}  ({time.time()-t0:.0f}s)')
        break
    a, u, p2, tr, av2, nz2 = got
    print(f'it{it}: a{a} via x_{u}  {cur[0]} -> {p2[0]}  nonzero {len(nz)}->{len(nz2)}',
          flush=True)
    v, cur, av, nz = tr, p2, av2, nz2
    if p2[0] > 39026:
        T.save(v, os.path.join(HERE, f'BREAK_{p2[0]}.json'))
        print(f'  *** BEATS THE DELIVERABLE: {p2[0]} saved', flush=True)
T.save(v, os.path.join(HERE, f'r4_{cur[0]}.json'))
print(f'FINAL {cur[0]} nonzero {nz}')
