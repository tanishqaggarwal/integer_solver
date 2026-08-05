"""Trace the EFFECTIVE mod-p computation of a variable: prune terms that vanish mod p
(products against the p-wire, products against currently-zero variables)."""
import pickle, collections, sys
import harness as H

P = 2**256 - 2**32 - 977
g = pickle.load(open('gates.pkl', 'rb'))
polys = pickle.load(open('polys.pkl', 'rb'))
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
definer = g['definer']
avars = [set(v for m in Pp for v in m) for Pp in polys]


def show(root, v, maxdepth=14, seen=None, depth=0, out=None):
    """Print the mod-p-surviving cone; mark leaves as FREE / p-wire / zero."""
    if out is None: out = []
    if seen is None: seen = set()
    pad = '  ' * depth
    a = definer.get(root)
    if a is None:
        tag = 'FREE'
        if v[root] % P == 0 and v[root] != 0: tag = 'FREE(=0 mod p)'
        out.append(f'{pad}x_{root} = {tag}  [val%p={v[root]%P}]')
        return out
    if root in seen or depth >= maxdepth:
        out.append(f'{pad}x_{root} ... (repeat/depth)')
        return out
    seen.add(root)
    Pp = polys[a]
    # which monomials survive mod p at the current point?
    live = []
    for m, c in Pp.items():
        if not m or (len(m) == 1 and m[0] == root):
            continue
        prod = c
        for u in m: prod *= v[u]
        live.append((m, c, prod % P))
    out.append(f'{pad}x_{root} <- [{a}] {src[a][:100]}   [val%p={v[root]%P}]')
    kids = set()
    for m, c, pv in live:
        for u in m:
            if u == root: continue
            kids.add(u)
    for u in sorted(kids):
        if v[u] % P == 0:
            out.append(f'{pad}  x_{u} = 0 mod p (dead branch)  [val={v[u]}]')
            continue
        show(u, v, maxdepth, seen, depth + 1, out)
    return out


if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    for r in [int(x) for x in sys.argv[1:]] or [12186]:
        print(f'======== effective mod-p cone of x_{r} (target K1%p = {v[2099]%P})')
        print('\n'.join(show(r, v)))
        print()
