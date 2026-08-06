"""S10 step 111: the repair engine, with the RIGHT potential function.

Score alone is a plateau: zeroing an atom often gains no equation because its
equations also contain another nonzero atom.  Rank states by
    (equations satisfied, -number of nonzero atoms, -total residue size)
so atom-count progress is never rejected.  Moves: zero-atom (2-level through a
gate's free input), and mod-p Newton on free inputs.
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

def pot(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    s = L.NEQ - len(L.failing_eqs(av))
    return (s, -len(nz), -sum(abs(av[a]).bit_length() for a in nz)), av, nz

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

def run(v, tag, iters=60, budget=2400):
    cur, av, nz = pot(v)
    print(f'{tag}: start {cur[0]} (nonzero {len(nz)})', flush=True)
    t0 = time.time()
    for it in range(iters):
        if time.time() - t0 > budget: break
        got = None
        for a in nz:
            for u, nv in moves(a, v, av):
                tr = list(v); tr[u] = nv
                ad.fwd(tr, rounds=6)
                p2, av2, nz2 = pot(tr)
                if p2 > cur: got = (a, u, p2, tr, av2, nz2); break
            if got: break
        if not got:
            print(f'  {tag} it{it}: stuck  score {cur[0]}  nonzero {nz}', flush=True)
            break
        a, u, p2, tr, av2, nz2 = got
        print(f'  {tag} it{it}: a{a} via x_{u}  score {cur[0]} -> {p2[0]}  '
              f'nonzero {len(nz)} -> {len(nz2)}', flush=True)
        v, cur, av, nz = tr, p2, av2, nz2
        if p2[0] > 39026:
            T.save(v, os.path.join(HERE, f'ENGINE_{p2[0]}.json'))
    T.save(v, os.path.join(HERE, f'engine_{tag}_{cur[0]}.json'))
    print(f'{tag} FINAL {cur[0]} nonzero {nz}', flush=True)
    return v, cur

if __name__ == '__main__':
    which = sys.argv[1]
    if which.startswith('br'):
        base = L.load(os.path.join(HERE, 'forward_state.json'))
        v = list(base); v[2081] = int(which[2]); v[4287] = int(which[3])
        ad.fwd(v, rounds=6)
    else:
        v = L.load(os.path.join(HERE, which))
    run(v, os.path.basename(which).replace('.json', ''))
