"""Derive gate candidates and an acyclic definer orientation over all atoms."""
import pickle, collections, sys

NVARS = 38748

def load():
    polys = pickle.load(open('polys.pkl','rb'))
    return polys

def gate_candidates(polys):
    """For each atom, the vars it can define: linear coeff +-1 and absent from all higher monomials."""
    cands = []   # atom -> list of vars
    for P in polys:
        lin = {}; higher = set()
        for m, c in P.items():
            if len(m) == 1: lin[m[0]] = c
            elif len(m) >= 2: higher.update(m)
        cands.append([t for t, c in lin.items() if c in (1, -1) and t not in higher])
    return cands

def orient(polys, cands, seed_free=(), forbid_def=()):
    """Greedy topological orientation: a gate defines its target once all other vars are ready."""
    var_gates = collections.defaultdict(list)
    for a, ts in enumerate(cands):
        for t in ts: var_gates[t].append(a)
    atom_vars = [set(v for m in P for v in m) for P in polys]
    ready = [False]*NVARS
    definer = {}
    forbid = set(forbid_def)
    # free inputs: no candidate definer at all
    for v in range(NVARS):
        if v not in var_gates: ready[v] = True
    for v in seed_free: ready[v] = True
    # count of unready vars per (atom,target)
    pending = collections.defaultdict(int)
    watch = collections.defaultdict(list)   # var -> list of (atom,target)
    q = collections.deque()
    for a, ts in enumerate(cands):
        for t in ts:
            if (a, t) in forbid: continue
            need = [u for u in atom_vars[a] if u != t and not ready[u]]
            pending[(a, t)] = len(need)
            for u in need: watch[u].append((a, t))
            if not need: q.append((a, t))
    order = []
    while q:
        a, t = q.popleft()
        if ready[t]: continue
        ready[t] = True; definer[t] = a; order.append(t)
        for key in watch[t]:
            pending[key] -= 1
            if pending[key] == 0: q.append(key)
    return definer, order, ready, var_gates

if __name__ == '__main__':
    polys = load()
    cands = gate_candidates(polys)
    ncand = sum(1 for c in cands if c)
    print(f'atoms with >=1 gate candidate: {ncand}/{len(polys)}')
    definer, order, ready, var_gates = orient(polys, cands)
    freeinp = [v for v in range(NVARS) if v not in var_gates]
    undef = [v for v in range(NVARS) if v in var_gates and v not in definer]
    print(f'free inputs (no candidate definer): {len(freeinp)}')
    print(f'defined by orientation: {len(definer)}   undefined-but-definable: {len(undef)}')
    used = set(definer.values())
    print(f'atoms used as definitions: {len(used)}   atoms left as CHECKS: {len(polys)-len(used)}')
    pickle.dump({'cands':cands,'definer':definer,'order':order,'freeinp':freeinp,'undef':undef},
                open('dag.pkl','wb'))
    for name, a in [('C1',22229),('C2',22231)]:
        t = [t for t,aa in definer.items() if aa==a]
        print(f'{name} atom {a}: candidates={cands[a]}, used-as-definer-of={t}')
