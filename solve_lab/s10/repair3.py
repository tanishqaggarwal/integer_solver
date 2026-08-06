"""S10 step 109: enriched repair INSIDE frame 2, with the seven p-checks
renormalised after every move (they are always satisfiable, so they cost nothing).
Remaining targets: the gadget cluster a7930 and a29539."""
import os, sys, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import DETACH, definer, ORDER, FREE, CHECKS, fwd, grad
P = ad.P
FORBID = {2081, 4287}
PINNED = {1329, 29854, 31864, 10903, 642, 28730, 7068, 9118, 8731, 17325, 9413}

def renorm(v):
    v[1329] = 5113045 * v[7075] * (v[9118] // P)
    v[29854] = 5113045 * v[7075] * v[9118]
    v[31864] = -v[7075] * v[8731]
    v[10903] = v[31864] // P
    v[642] = v[17325] * P
    v[28730] = v[9413] * P
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    return v

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def moves(a, v, av):
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID or w in PINNED: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in FREE:
            out.append((w, tgt))
        else:
            d = definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE or u in FORBID or u in PINNED: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: out.append((u, nv))
    vm = [x % P for x in v]
    r = av[a] % P
    if r:
        try: g = grad(a, vm)
        except Exception: g = {}
        cand = sorted(((len(L.var_atoms[u]), u, d) for u, d in g.items()
                       if u not in FORBID and u not in PINNED and d % P))
        for _, u, d in cand[:16]:
            out.append((u, v[u] + (-r * pow(d, -1, P)) % P))
    return out

v = L.load(os.path.join(HERE, 'construct_39004.json')) if os.path.exists(
    os.path.join(HERE, 'construct_39004.json')) else None
if v is None:
    from construct2 import build, base
    v = renorm(build(base))
best = score(v)
print(f'frame-2 constructed state: score {best}', flush=True)
t0 = time.time()
for it in range(30):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    got = None
    for a in nz:
        for u, nv in moves(a, v, av):
            tr = list(v); tr[u] = nv
            renorm(tr)
            s = score(tr)
            if s > best: got = (a, u, s, tr); break
        if got: break
    if not got:
        print(f'it{it}: stuck at {best}; nonzero {nz}  ({time.time()-t0:.0f}s)')
        break
    a, u, s, tr = got
    print(f'it{it}: a{a} via x_{u}  {best} -> {s}', flush=True)
    v, best = tr, s
    T.save(v, os.path.join(HERE, f'f2_{s}.json'))
av = L.all_atom_values(v)
print(f'FINAL score {best}  nonzero {[a for a in range(L.NA) if av[a]]}')
