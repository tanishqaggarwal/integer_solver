"""Target 39,024: close chain 2 exactly, leave chain 1 (atom 22229) broken.
Relax 22229 from the gate set so x_7068 is a FREE input and C1 is NOT a constraint."""
import pickle, collections, heapq, sys, time, os
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0, S9); os.chdir(S9)
import harness as H
P = 2**256-2**32-977
NV = 38748
g = pickle.load(open('gates.pkl','rb'))
polys = pickle.load(open('polys.pkl','rb'))
topo = pickle.load(open('topo.pkl','rb'))['topo']
roots = pickle.load(open('roots.pkl','rb'))
definer = dict(g['definer']); atom_out = dict(g['atom_out'])
RELAX = 22229                      # <-- chain 1 left free
del atom_out[RELAX]; del definer[7068]
avars = [set(v for m in Pp for v in m) for Pp in polys]
var_atoms = collections.defaultdict(list)
for a, s in enumerate(avars):
    for u in s: var_atoms[u].append(a)
tidx = {v: i for i, v in enumerate(topo)}; BIG = len(topo)+10
checks = [a for a in range(len(polys)) if a not in atom_out and a != RELAX]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}
freeinp = [x for x in range(NV) if x not in definer]

def evalpoly(Pp, v):
    s = 0
    for m, c in Pp.items():
        t = c
        for u in m: t *= v[u]
        s += t
    return s

def solve_for(a, t, v):
    Pp = polys[a]; c = 0
    for m, cc in Pp.items():
        if len(m) == 1 and m[0] == t: c += cc
        elif t in m: return None
    if c == 0: return None
    old = v[t]; v[t] = 0; rest = evalpoly(Pp, v); v[t] = old
    if rest % c: return None
    return -rest // c

def ripple(v, seeds):
    changed = {}; pq = []
    for u, val in seeds.items():
        if v[u] == val: continue
        v[u] = val; changed[u] = val; heapq.heappush(pq, (tidx.get(u, BIG), u))
    while pq:
        _, u = heapq.heappop(pq)
        for a in var_atoms[u]:
            if a not in atom_out: continue
            c, t = atom_out[a]
            if t == u: continue
            nv = solve_for(a, t, v)
            if nv is None or nv == v[t]: continue
            v[t] = nv; changed[t] = nv; heapq.heappush(pq, (tidx.get(t, BIG), t))
    return changed

def resid(v): return {a: evalpoly(Pp, v) for a, Pp in rp.items()}

if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    print('baseline atom22229 =', evalpoly(polys[22229], v))
    ripple(v, {4432: v[19964] + v[28730]})          # close C2 exactly
    print('after closing C2: atom22231 =', evalpoly(polys[22231], v))
    base = resid(v)
    nz = {a: x for a, x in base.items() if x}
    print('nonzero checks (chain-1 relaxed):', sorted(nz))
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'EQUATIONS now: {len(codes)-len(f)}/{len(codes)} ({len(f)} failing) {f}')
    pickle.dump({'checks':checks}, open('chain2/meta.pkl','wb'))
    H.save_assignment(v, 'chain2/state_c2.json')
