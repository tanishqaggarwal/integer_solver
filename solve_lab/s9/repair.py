"""Iterative ripple-repair: fix broken atoms via free inputs, re-ripple, repeat."""
import pickle, collections, heapq, json, sys
import harness as H
NVARS = 38748

g = pickle.load(open('gates.pkl','rb'))
polys = pickle.load(open('polys.pkl','rb'))
topo = pickle.load(open('topo.pkl','rb'))['topo']
d = pickle.load(open('atoms.pkl','rb')); src = d['atom_src']
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

def solve_for(a, t, v):
    """Value of var t making atom a vanish (t must appear linearly), or None."""
    P = polys[a]
    c = 0; 
    for m, cc in P.items():
        if len(m) == 1 and m[0] == t: c += cc
        elif t in m: return None          # nonlinear in t
    if c == 0: return None
    old = v[t]; v[t] = 0
    rest = evalpoly(P, v)
    v[t] = old
    if rest % c: return None
    return -rest // c

def ripple(v, seeds):
    changed = {}; pq = []
    for u, val in seeds.items():
        if v[u] == val: continue
        v[u] = val; changed[u] = val
        heapq.heappush(pq, (tidx.get(u, BIG), u))
    steps = 0
    while pq:
        _, u = heapq.heappop(pq); steps += 1
        if steps > 500000: break
        for a in var_atoms[u]:
            if a not in atom_out: continue
            c, t = atom_out[a]
            if t == u: continue
            nv = solve_for(a, t, v)
            if nv is None or nv == v[t]: continue
            v[t] = nv; changed[t] = nv
            heapq.heappush(pq, (tidx.get(t, BIG), t))
    return changed, steps

def nonzero_atoms(v):
    return [a for a, P in enumerate(polys) if evalpoly(P, v) != 0]

def repair_loop(v, rounds=40, verbose=True):
    freeset = set(x for x in range(NVARS) if x not in definer)
    history = []
    for r in range(rounds):
        nz = nonzero_atoms(v)
        history.append(len(nz))
        if verbose: print(f'--- round {r}: nonzero atoms = {len(nz)}  {nz[:12]}')
        if not nz: return True, history
        seeds = {}
        stuck = []
        for a in nz:
            best = None
            for t in sorted(avars[a]):
                if t not in freeset: continue
                nv = solve_for(a, t, v)
                if nv is None: continue
                # prefer the free var with the fewest atom occurrences
                k = len(var_atoms[t])
                if best is None or k < best[0]: best = (k, t, nv)
            if best is None: stuck.append(a)
            else:
                k, t, nv = best
                if t in seeds and seeds[t] != nv: stuck.append(a)
                else: seeds[t] = nv
        if verbose: print(f'    repairs: {len(seeds)} free vars, stuck atoms: {len(stuck)} {stuck[:8]}')
        if not seeds: return False, history
        changed, steps = ripple(v, seeds)
        if verbose: print(f'    ripple changed {len(changed)} vars')
    return False, history

if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    seeds = {7068: v[2099] + 7376877*v[642], 4432: v[19964] + v[28730]}
    ripple(v, seeds)
    ok, hist = repair_loop(v)
    print('converged:', ok, 'history:', hist)
    codes, _ = H.load_equations()
    fails = H.evaluate(codes, v)
    print(f'EQUATIONS satisfied: {len(codes)-len(fails)}/{len(codes)} fails={fails[:20]}')
    H.save_assignment(v, 'repair_out.json')
