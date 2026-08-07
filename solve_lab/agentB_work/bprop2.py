"""Agent B: directed repair.  Change a set of variables, then propagate the change
through the gate graph BFS-style: each violated gate is repaired by assigning ONE
not-yet-settled variable (preferring sink-like wires).  A variable is settled once
assigned, so the propagation is a DAG walk and cannot oscillate.
"""
import pickle, collections, json, sys, time

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
NV = 38748

occ = collections.defaultdict(list)
fvars = []
for i, p in enumerate(facs):
    vs = set()
    for m in p: vs.update(m)
    fvars.append(vs)
    for v in vs: occ[v].append(i)
cands = []
for p in facs:
    sq = set()
    for m in p:
        if len(m) == 2 and m[0] == m[1]: sq.add(m[0])
    allv = set()
    for m in p: allv.update(m)
    cands.append(sorted(allv - sq))
atom_of_fac = collections.defaultdict(list)
for i, a in enumerate(atoms):
    for f in a: atom_of_fac[f].append(i)

def load(path):
    val = [0]*NV
    for k, v in json.load(open(path)).items(): val[int(k.split('_')[1])] = int(v)
    return val

def fval(f, val):
    t = 0
    for m, c in facs[f].items():
        z = c
        for v in m:
            z *= val[v]
            if z == 0: break
        t += z
    return t

def solve_for(f, v, val):
    """solve gate f for v (v must occur to degree <=1 in every monomial)."""
    k = 0; a = 0
    for m, c in facs[f].items():
        if v in m:
            z = c
            for u in m:
                if u == v: continue
                z *= val[u]
                if z == 0: break
            a += z
        else:
            z = c
            for u in m:
                z *= val[u]
                if z == 0: break
            k += z
    if a == 0 or (-k) % a: return None
    return (-k)//a

def score(val):
    fv = [fval(f, val) for f in range(len(facs))]
    ok = 0; fail = []
    for e, (sc, L, k) in enumerate(eqs):
        s = 0
        for c, a in L:
            t = c
            for f in atoms[a]:
                t *= fv[f]
                if t == 0: break
            s += t
        if s == 0: ok += 1
        else: fail.append(e)
    return ok, fail, fv

def propagate(val, forced, base_fv, maxsteps=200000, verbose=False):
    """forced: dict var->new value.  base_fv: gate values BEFORE the change
       (a gate that was already nonzero before stays 'allowed')."""
    allowed = set(f for f in range(len(facs)) if base_fv[f] != 0)
    settled = set(forced)
    for v, x in forced.items(): val[v] = x
    q = collections.deque()
    seen = set()
    for v in forced:
        for f in occ[v]:
            if f not in seen: q.append(f); seen.add(f)
    steps = 0; unfixable = []
    while q:
        f = q.popleft(); seen.discard(f); steps += 1
        if steps > maxsteps: break
        if fval(f, val) == 0: continue
        if f in allowed: continue
        cs = [v for v in cands[f] if v not in settled]
        if not cs:
            unfixable.append(f); continue
        # prefer sink-like wires: fewest occurrences
        cs.sort(key=lambda v: (len(occ[v]), v))
        done = False
        for v in cs:
            x = solve_for(f, v, val)
            if x is None: continue
            val[v] = x; settled.add(v); done = True
            for g in occ[v]:
                if g != f and g not in seen: q.append(g); seen.add(g)
            break
        if not done: unfixable.append(f)
    return val, unfixable, settled

if __name__ == '__main__':
    base = sys.argv[1]
    forced = {}
    for kv in sys.argv[2].split(','):
        k, v = kv.split('='); forced[int(k)] = int(v)
    val = load(base)
    ok0, fail0, fv0 = score(val)
    print("base %d  failing=%s" % (ok0, fail0))
    val, unf, settled = propagate(val, forced, fv0)
    ok, fail, fv = score(val)
    print("after propagation: %d  settled=%d  unfixable=%d  failing=%s" % (ok, len(settled), len(unf), fail[:15]))
    if len(sys.argv) > 3:
        json.dump({('x_%d'%i): val[i] for i in range(NV)}, open(sys.argv[3], 'w'))
        print("wrote", sys.argv[3])
