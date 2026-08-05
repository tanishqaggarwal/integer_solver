"""Topological forward evaluation from free inputs using the canonical gate orientation."""
import pickle, collections, sys, json
import harness as H
NVARS = 38748

def build():
    g = pickle.load(open('gates.pkl','rb'))
    polys = pickle.load(open('polys.pkl','rb'))
    definer, atom_out, free = g['definer'], g['atom_out'], g['free']
    avars = [set(v for m in P for v in m) for P in polys]
    # topological sort over defined vars
    indeg = {}; users = collections.defaultdict(list)
    for v, a in definer.items():
        deps = [u for u in avars[a] if u != v and u in definer]
        indeg[v] = len(deps)
        for u in deps: users[u].append(v)
    q = collections.deque([v for v in definer if indeg[v] == 0])
    topo = []
    while q:
        v = q.popleft(); topo.append(v)
        for w in users[v]:
            indeg[w] -= 1
            if indeg[w] == 0: q.append(w)
    return polys, definer, atom_out, free, topo, avars

def evalpoly(P, v):
    s = 0
    for m, c in P.items():
        t = c
        for u in m: t *= v[u]
        s += t
    return s

if __name__ == '__main__':
    polys, definer, atom_out, free, topo, avars = build()
    print(f'topological order covers {len(topo)}/{len(definer)} defined vars '
          f'({len(definer)-len(topo)} in cycles)')
    pickle.dump({'topo':topo}, open('topo.pkl','wb'))
