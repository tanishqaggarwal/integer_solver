"""Expand a boolean OR/AND tree down to non-boolean leaves."""
import sys, collections
import dlib as L

P = L.P

def _pair(u):
    """If x_u = a + b, return (a,b)."""
    a = L.definer.get(u)
    if a is None: return None
    p = L.polys[a]
    if p.get((u,)) != 1: return None
    rest = {m: -c for m, c in p.items() if m != (u,)}
    if len(rest) == 2 and all(len(m) == 1 and c == 1 for m, c in rest.items()):
        return tuple(m[0] for m in rest)
    return None


def _prod(u):
    """If x_u = a*b, return (a,b)."""
    a = L.definer.get(u)
    if a is None: return None
    p = L.polys[a]
    if p.get((u,)) != 1: return None
    rest = {m: -c for m, c in p.items() if m != (u,)}
    if len(rest) == 1:
        m, c = list(rest.items())[0]
        if len(m) == 2 and c == 1:
            return m
    return None


def kind(u):
    """Classify definition of u: ('copy',w) ('or',a,b) ('and',a,b) ('sum',...) ('free',) ('other',)"""
    a = L.definer.get(u)
    if a is None:
        return ('free',)
    p = L.polys[a]
    # u - w
    terms = {m: c for m, c in p.items()}
    # normalise so coeff of (u,) is 1
    cu = terms.get((u,), 0)
    if cu == 0:
        return ('other', a)
    rest = {m: -c for m, c in terms.items() if m != (u,)}
    if cu != 1:
        return ('other', a)
    if len(rest) == 1 and list(rest)[0] != () and len(list(rest)[0]) == 1:
        return ('copy', list(rest)[0][0])
    if len(rest) == 1 and len(list(rest)[0]) == 2:
        m = list(rest)[0]
        if rest[m] == 1:
            return ('and', m[0], m[1])
    # u - (A - B) with A = a+b, B = a*b  -> OR(a,b)
    if len(rest) == 2:
        ks = sorted(rest.items(), key=lambda t: -t[1])
        if len(ks[0][0]) == 1 and len(ks[1][0]) == 1 and ks[0][1] == 1 and ks[1][1] == -1:
            A, B = ks[0][0][0], ks[1][0][0]
            ka, kb = _pair(A), _prod(B)
            if ka and kb and set(ka) == set(kb):
                return ('or', ka[0], ka[1])
    if len(rest) == 3:
        # a + b - a*b  (OR)
        lin = [m for m in rest if len(m) == 1]
        quad = [m for m in rest if len(m) == 2]
        if len(lin) == 2 and len(quad) == 1 and rest[quad[0]] == -1 and all(rest[m] == 1 for m in lin):
            if set(quad[0]) == set(x[0] for x in lin):
                return ('or', lin[0][0], lin[1][0])
    return ('other', a)


def expand(u, v, ind=0, seen=None, out=None):
    if seen is None: seen = set()
    if out is None: out = []
    k = kind(u)
    if k[0] == 'copy':
        return expand(k[1], v, ind, seen, out)
    if k[0] in ('or', 'and'):
        print('  '*ind + f'x_{u} = {k[0].upper()}(x_{k[1]}, x_{k[2]})  val={v[u]}')
        expand(k[1], v, ind+1, seen, out)
        expand(k[2], v, ind+1, seen, out)
        return out
    a = L.definer.get(u)
    src = 'FREE' if a is None else L.atom_src[a][:150]
    print('  '*ind + f'LEAF x_{u} val={v[u]} :: {src}')
    out.append(u)
    return out


if __name__ == '__main__':
    v = L.load(sys.argv[1])
    for u in [int(x) for x in sys.argv[2:]]:
        print('='*100)
        leaves = expand(u, v)
        print('leaves:', leaves)
