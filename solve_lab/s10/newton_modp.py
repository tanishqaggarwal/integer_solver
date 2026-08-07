"""S10 step 95: mod-p NEWTON moves on free inputs + p-handle absorption.

For a residual atom a whose only exact repair is a p-quantised handle h, the
requirement is not a value but  a == 0 (mod p).  Take the exact AD derivative
d = da/du (mod p) for a free input u, solve  a + d*delta == 0 (mod p)  for delta,
and move u by delta.  Then the handle absorbs a/p exactly over Z.

No search so far could propose this: shifting u by delta zeroes NO atom on its
own, so every zero-this-atom move generator is blind to it.
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

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def newton_moves(a, v, av, topn=40):
    """(u, delta) shifting free input u so that atom a becomes 0 mod p."""
    vm = [x % P for x in v]
    g = ad.grad(a, vm)
    r = av[a] % P
    if r == 0: return []
    out = []
    for u, d in g.items():
        if u in FORBID or d % P == 0: continue
        delta = (-r * pow(d, -1, P)) % P
        out.append((u, delta, len(L.var_atoms[u])))
    out.sort(key=lambda t: t[2])          # fewest consumers first
    return [(u, d) for u, d, _ in out[:topn]]

if __name__ == '__main__':
    W = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v0 = L.load(W)
    s0 = score(v0)
    av0 = L.all_atom_values(v0)
    nz0 = [a for a in range(L.NA) if av0[a]]
    print(f'delivered witness: score {s0}, nonzero {nz0}')
    for a in nz0:
        mvs = newton_moves(a, v0, av0)
        print(f'\na{a} ({len(L.atom2eq[a])} eqs) residue mod p nonzero; '
              f'{len(mvs)} single-variable mod-p Newton moves available')
        best = None
        for u, delta in mvs[:14]:
            v = list(v0); v[u] = v[u] + delta
            ad.fwd(v, rounds=6)
            av = L.all_atom_values(v)
            s = score(v)
            nz = [b for b in range(L.NA) if av[b]]
            hit = av[a] % P == 0
            if best is None or s > best[0]: best = (s, u, delta, nz, hit)
            print(f'   x_{u:<6} (consumers {len(L.var_atoms[u]):>2}) -> score {s:>6} '
                  f'({s-s0:+d})  a{a} mod p zero? {hit}  nonzero {nz[:9]}')
        if best: print(f'   BEST {best[0]} via x_{best[1]}')
