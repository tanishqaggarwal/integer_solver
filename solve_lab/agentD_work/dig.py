import sys, collections
import dlib as L
import fwd

P = L.P

def d1(u, v):
    a = L.definer.get(u)
    if a is None:
        return f'x_{u} = FREE'
    return f'x_{u} :: a{a} = {L.atom_src[a]}'

def tree(u, v, depth, seen=None, ind=0, maxlen=200):
    if seen is None: seen = set()
    val = v[u]
    sval = str(val)
    if len(sval) > 30: sval = sval[:12]+'...'+f'({len(sval)}d, mod p={val%P if val else 0})'[:60]
    print('  '*ind + d1(u, v)[:maxlen] + '   |val=' + sval)
    if depth <= 0 or u in seen: return
    seen.add(u)
    a = L.definer.get(u)
    if a is None: return
    for w in sorted(L.avars[a]):
        if w != u:
            tree(w, v, depth-1, seen, ind+1, maxlen)

if __name__ == '__main__':
    v = L.load(sys.argv[1]) if len(sys.argv) > 1 else None
    if v is None:
        v = L.load('D_state1.json')
    args = [int(a) for a in sys.argv[2:]] or [15298]
    depth = 4
    for u in args:
        print('='*100)
        tree(u, v, depth)
