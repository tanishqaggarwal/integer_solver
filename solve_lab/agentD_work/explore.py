import sys, collections
import dlib as L

P = L.P

def isfree(u):
    return u in L.freeset

def desc(u, v=None):
    d = L.definer.get(u)
    s = f'x_{u}: '
    if d is None:
        s += 'FREE'
    else:
        s += f'def a{d} = {L.atom_src[d][:140]}'
    if v is not None:
        s += f'   val={v[u]}'
    return s

def tree(u, v, depth=2, seen=None, ind=0):
    if seen is None: seen = set()
    print('  '*ind + desc(u, v))
    if depth == 0 or u in seen: return
    seen.add(u)
    d = L.definer.get(u)
    if d is None: return
    for w in sorted(L.avars[d]):
        if w != u:
            tree(w, v, depth-1, seen, ind+1)

if __name__ == '__main__':
    import fwd
    v = L.load('../best/new_instance_partial_39026.json')
    for t in L.definer: v[t] = 0
    fwd.forward(v)
    for u in [9118, 8731, 24548, 25442, 14853, 1308, 7075, 7927, 29967, 29854, 31864]:
        print(desc(u, v))
        print('    is free:', isfree(u), ' value mod p:', v[u] % P if v[u] else 0)
    print()
    print('--- a7930 tree ---')
    tree(24548, v, 3)
    print()
    tree(25442, v, 3)
