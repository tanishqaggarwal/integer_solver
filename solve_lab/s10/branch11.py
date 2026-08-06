"""S10 step 110: full enriched engine in branch (1,1).

With x_7075 = 0 the two congruences  p | x_9118  and  p | x_8731  VANISH -- a35759
and a35761 become x_29854 = 0 and x_31864 = 0, satisfiable through the solo
handles x_1329 and x_10903 at zero cost.  That removes two of the constraints
fighting the gadget cluster.  Price: the (1,1) branch activates x_21279 = 1,
switching on a19088, a22233, a22235.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
FORBID = {2081, 4287}

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def moves(a, v, av, nnewton=20):
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in FREE: out.append((w, tgt))
        else:
            d = definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE or u in FORBID: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: out.append((u, nv))
    r = av[a] % P
    if r:
        vm = [x % P for x in v]
        try: g = ad.grad(a, vm)
        except Exception: g = {}
        cand = sorted((len(L.var_atoms[u]), u, d) for u, d in g.items()
                      if u not in FORBID and d % P)
        for _, u, d in cand[:nnewton]:
            out.append((u, v[u] + (-r * pow(d, -1, P)) % P))
    return out

base = L.load(os.path.join(HERE, 'forward_state.json'))
v = list(base); v[2081] = 1; v[4287] = 1
ad.fwd(v, rounds=6)
best = score(v)
print(f'branch(1,1) raw score {best}', flush=True)
t0 = time.time()
for it in range(40):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    got = None
    for a in nz:
        for u, nv in moves(a, v, av):
            tr = list(v); tr[u] = nv
            ad.fwd(tr, rounds=6)
            s = score(tr)
            if s > best: got = (a, u, s, tr); break
        if got: break
    if not got:
        print(f'it{it}: stuck at {best}; nonzero {nz}  ({time.time()-t0:.0f}s)')
        break
    a, u, s, tr = got
    print(f'it{it}: a{a} via x_{u}  {best} -> {s}  ({time.time()-t0:.0f}s)', flush=True)
    v, best = tr, s
    T.save(v, os.path.join(HERE, f'b11_{s}.json'))
av = L.all_atom_values(v)
print(f'FINAL {best}  nonzero {[a for a in range(L.NA) if av[a]]}')
