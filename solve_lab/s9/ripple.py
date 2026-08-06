"""Set C1/C2 gate outputs to canonical values and forward-ripple through the gate DAG."""
import pickle, collections, heapq, json, sys
import harness as H
NVARS = 38748

g = pickle.load(open('gates.pkl','rb'))
polys = pickle.load(open('polys.pkl','rb'))
topo = pickle.load(open('topo.pkl','rb'))['topo']
definer, atom_out = g['definer'], g['atom_out']
avars = [set(v for m in P for v in m) for P in polys]
var_atoms = collections.defaultdict(list)
for a, s in enumerate(avars):
    for u in s: var_atoms[u].append(a)
tidx = {v: i for i, v in enumerate(topo)}
BIG = len(topo) + 10

def evalpoly(P, v):
    s = 0
    for m, c in P.items():
        t = c
        for u in m: t *= v[u]
        s += t
    return s

def solve_out(a, v):
    """Return the value of atom a's output var that makes atom a vanish, or None."""
    c, t = atom_out[a]
    old = v[t]; v[t] = 0
    rest = evalpoly(polys[a], v)
    v[t] = old
    if rest % c: return None
    return -rest // c

def ripple(v, seeds, maxsteps=2_000_000, verbose=True):
    """seeds: {var: newvalue}. Returns (changed_vars, unresolved_atoms)."""
    changed = {}
    pq = []
    for u, val in seeds.items():
        v[u] = val; changed[u] = val
        heapq.heappush(pq, (tidx.get(u, BIG), u))
    seen = set()
    bad = []
    steps = 0
    while pq and steps < maxsteps:
        _, u = heapq.heappop(pq); steps += 1
        for a in var_atoms[u]:
            if a not in atom_out: continue
            c, t = atom_out[a]
            if t == u: continue
            nv = solve_out(a, v)
            if nv is None:
                bad.append(('nodiv', a)); continue
            if nv != v[t]:
                v[t] = nv; changed[t] = nv
                heapq.heappush(pq, (tidx.get(t, BIG), t))
    return changed, bad, steps

if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    seeds = {}
    # canonical: x_7068 := x_2099 + 7376877*x_642 ; x_4432 := x_19964 + x_28730
    seeds[7068] = v[2099] + 7376877*v[642]
    seeds[4432] = v[19964] + v[28730]
    print('delta x_7068 =', seeds[7068]-v[7068])
    print('delta x_4432 =', seeds[4432]-v[4432])
    changed, bad, steps = ripple(v, seeds)
    print(f'ripple: {len(changed)} vars changed, {len(bad)} non-divisible gates, {steps} steps')
    nz = [a for a, P in enumerate(polys) if evalpoly(P, v) != 0]
    print(f'nonzero atoms after ripple: {len(nz)}  {nz[:30]}')
    codes, _ = H.load_equations()
    fails = H.evaluate(codes, v)
    print(f'EQUATIONS satisfied: {len(codes)-len(fails)}/{len(codes)}  fails={fails[:30]}')
    if len(fails) < 11:
        H.save_assignment(v, 'ripple_out.json'); print('saved ripple_out.json')
