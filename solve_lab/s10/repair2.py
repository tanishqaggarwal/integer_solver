"""S10 step 93: two-level repair engine, controls FORBIDDEN.

set_var(w, target): realise x_w = target either directly (w free) or through a
free input of w's definer.  This is exactly the handle-based repair lib.ripple
cannot see.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
FORBID = {2081, 4287}

def solve_for_target(a, w, v, target):
    """value of free u in atom a making atom a vanish when x_w = target."""
    vv = list(v); vv[w] = target
    out = []
    for u in sorted(set(L.avars[a])):
        if u == w or u not in FREE or u in FORBID: continue
        nv = T.solve_lin(a, u, vv)
        if nv is not None: out.append((u, nv))
    return out

def candidates(a, v):
    """(var, newvalue) moves that would zero atom a."""
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        nv = T.solve_lin(a, w, v)
        if nv is None or nv == v[w]: continue
        if w in FREE:
            out.append(((w, nv),))
        else:
            d = definer.get(w)
            if d is None: continue
            for u, uv in solve_for_target(d, w, v, nv):
                out.append(((u, uv),))
    return out

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def run(v, tag, iters=25):
    best = score(v)
    print(f'{tag}: start score {best}')
    for it in range(iters):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        moved = None
        for a in nz:
            for mv in candidates(a, v):
                trial = list(v)
                for u, uv in mv: trial[u] = uv
                ad.fwd(trial, rounds=6)
                s = score(trial)
                if s > best:
                    moved = (a, mv, s, trial); break
            if moved: break
        if not moved:
            print(f'  it{it}: stuck at {best}; nonzero {nz}')
            return v, best
        a, mv, s, trial = moved
        print(f'  it{it}: a{a} via {[f"x_{u}" for u, _ in mv]}  {best} -> {s}')
        v, best = trial, s
    return v, best

if __name__ == '__main__':
    base = L.load(os.path.join(HERE, 'forward_state.json'))
    for b1, b2 in [(1, 1), (1, 0)]:
        v = list(base); v[2081] = b1; v[4287] = b2
        ad.fwd(v, rounds=6)
        v, s = run(v, f'branch({b1},{b2})')
        T.save(v, os.path.join(HERE, f'br{b1}{b2}.json'))
        print(f'branch({b1},{b2}) FINAL score {s}\n')
