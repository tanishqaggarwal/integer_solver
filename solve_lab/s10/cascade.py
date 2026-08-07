"""S11 step 3: the REPAIR CASCADE.

Every search so far hill-climbed: a move had to improve the score immediately.
But this instance is a chain -- fixing gadget G moves free input u, which breaks
every gadget downstream of u, and only after the whole chain is repaired does the
score go up.  A hill-climber can never walk that.

Cascade instead: apply the seed move, then repeatedly repair EVERY nonzero atom,
preferring repair variables with the fewest other consumers (a solo handle has
none, so it terminates the chain).  Accept regardless of score; measure only at
the end.
"""
import os, sys, time, collections
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

def repairs(a, v):
    """(cost, var, newvalue) moves that zero atom a; cost = other consumers."""
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in FREE:
            out.append((len(L.var_atoms[w]) - 1, w, tgt))
        else:
            d = definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE or u in FORBID: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None:
                    out.append((len(L.var_atoms[u]) - 1, u, nv))
    out.sort()
    return out

def newton_repairs(a, v, av, topn=6):
    r = av[a] % P
    if not r: return []
    vm = [x % P for x in v]
    try: g = ad.grad(a, vm)
    except Exception: return []
    c = sorted((len(L.var_atoms[u]) - 1, u, d) for u, d in g.items()
               if u not in FORBID and d % P)
    return [(k, u, v[u] + (-r * pow(d, -1, P)) % P) for k, u, d in c[:topn]]

def cascade(v, seed, maxsteps=40, verbose=False):
    v = list(v)
    u0, nv0 = seed
    v[u0] = nv0
    ad.fwd(v, rounds=6)
    seen = set()
    best = (score(v), list(v))
    for step in range(maxsteps):
        av = L.all_atom_values(v)
        nz = tuple(a for a in range(L.NA) if av[a])
        if not nz:
            return (L.NEQ, v, step)
        if nz in seen: break
        seen.add(nz)
        # repair the atom with the cheapest available terminating move
        cands = []
        for a in nz:
            for k, u, x in repairs(a, v)[:3]:
                cands.append((k, a, u, x))
            for k, u, x in newton_repairs(a, v, av):
                cands.append((k + 100, a, u, x))
        if not cands: break
        cands.sort()
        k, a, u, x = cands[0]
        v[u] = x
        ad.fwd(v, rounds=6)
        s = score(v)
        if s > best[0]: best = (s, list(v))
        if verbose:
            print(f'    step{step}: a{a} via x_{u} (cons {k}) -> {s} '
                  f'nonzero {len(nz)}', flush=True)
    return (best[0], best[1], -1)

if __name__ == '__main__':
    v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
    av0 = L.all_atom_values(v0)
    print(f'seed state {score(v0)}  nonzero {[a for a in range(L.NA) if av0[a]]}',
          flush=True)
    t0 = time.time()
    best = (score(v0), None, None)
    for a in (21617, 29539):
        for k, u, x in (repairs(a, v0)[:6] + newton_repairs(a, v0, av0, topn=10)):
            s, w, st = cascade(v0, (u, x))
            tag = 'SOLVED' if st >= 0 and s == L.NEQ else ''
            print(f'  seed a{a} via x_{u} (cons {k}) -> {s} {tag}', flush=True)
            if s > best[0]:
                best = (s, u, a)
                T.save(w, os.path.join(HERE, f'casc_{s}.json'))
                print(f'    *** NEW BEST {s} saved', flush=True)
    print(f'\nBEST {best[0]} via x_{best[1]} on a{best[2]}  ({time.time()-t0:.0f}s)')
