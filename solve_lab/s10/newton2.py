"""S10 step 96: mod-p Newton on the canonical frame's two surviving gadgets."""
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

def close_handles(v):
    """after fwd, try to zero every nonzero CHECK through a free var of a gate."""
    for _ in range(6):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        did = False
        for a in nz:
            for w in sorted(set(L.avars[a])):
                if w in FORBID: continue
                tgt = T.solve_lin(a, w, v)
                if tgt is None or tgt == v[w]: continue
                cands = []
                if w in FREE:
                    cands.append((w, tgt))
                else:
                    d = definer.get(w)
                    if d is not None:
                        vv = list(v); vv[w] = tgt
                        for u in sorted(set(L.avars[d])):
                            if u == w or u not in FREE or u in FORBID: continue
                            nv = T.solve_lin(d, u, vv)
                            if nv is not None: cands.append((u, nv))
                for u, nv in cands:
                    tr = list(v); tr[u] = nv
                    ad.fwd(tr, rounds=6)
                    if score(tr) > score(v):
                        v = tr; did = True; break
                if did: break
            if did: break
        if not did: break
    return v

v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
s0 = score(v0)
av0 = L.all_atom_values(v0)
print(f'canonical frame: score {s0}, nonzero {[a for a in range(L.NA) if av0[a]]}')
best = (s0, None, v0)
for a in [21617, 29539]:
    mvs = newton_moves(a, v0, av0, topn=60)
    print(f'\na{a}: {len(mvs)} mod-p Newton moves')
    for u, delta in mvs[:18]:
        v = list(v0); v[u] = v[u] + delta
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        hit = av[a] % P == 0
        v = close_handles(v)
        s = score(v)
        av = L.all_atom_values(v)
        nz = [b for b in range(L.NA) if av[b]]
        print(f'   x_{u:<6} (cons {len(L.var_atoms[u]):>2}) -> score {s:>6} ({s-s0:+d})'
              f'  a{a}=0 mod p? {hit}  nonzero {nz[:10]}')
        if s > best[0]:
            best = (s, (a, u), v)
            T.save(v, os.path.join(HERE, f'newton_{s}.json'))
            print(f'      *** NEW BEST {s} saved')
print(f'\nBEST {best[0]} via {best[1]}')
