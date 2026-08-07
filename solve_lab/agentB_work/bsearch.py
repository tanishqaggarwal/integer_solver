"""Agent B: randomised directed propagation with restarts + incremental scoring."""
import pickle, collections, json, sys, time, random

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
fac_atoms = collections.defaultdict(list)
for i, a in enumerate(atoms):
    for f in set(a): fac_atoms[f].append(i)
atom_eqs = collections.defaultdict(list)
for e, (sc, L, k) in enumerate(eqs):
    for c, a in L: atom_eqs[a].append(e)

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

def full_state(val):
    fv = [fval(f, val) for f in range(len(facs))]
    av = []
    for a in atoms:
        t = 1
        for f in a:
            t *= fv[f]
            if t == 0: break
        av.append(t)
    eqv = []
    for sc, L, k in eqs:
        eqv.append(sum(c*av[a] for c, a in L))
    return fv, av, eqv

def incr_score(val, changed, fv, av, eqv):
    """recompute only what `changed` (set of vars) touches; returns new (fv,av,eqv,ok)"""
    fv = list(fv); av = list(av); eqv = list(eqv)
    dirty_f = set()
    for v in changed: dirty_f.update(occ[v])
    dirty_a = set()
    for f in dirty_f:
        fv[f] = fval(f, val); dirty_a.update(fac_atoms[f])
    dirty_e = set()
    for a in dirty_a:
        t = 1
        for f in atoms[a]:
            t *= fv[f]
            if t == 0: break
        av[a] = t; dirty_e.update(atom_eqs[a])
    for e in dirty_e:
        eqv[e] = sum(c*av[a] for c, a in eqs[e][1])
    ok = sum(1 for x in eqv if x == 0)
    return fv, av, eqv, ok

def propagate(val, forced, allowed, rnd, maxsteps=100000):
    settled = set(forced)
    changed = set(forced)
    for v, x in forced.items(): val[v] = x
    q = collections.deque(); seen = set()
    for v in forced:
        for f in occ[v]:
            if f not in seen: q.append(f); seen.add(f)
    unf = []
    steps = 0
    while q:
        f = q.popleft(); seen.discard(f); steps += 1
        if steps > maxsteps: break
        if f in allowed: continue
        if fval(f, val) == 0: continue
        cs = [v for v in cands[f] if v not in settled]
        rnd.shuffle(cs)
        done = False
        for v in cs:
            x = solve_for(f, v, val)
            if x is None or x == val[v]: continue
            val[v] = x; settled.add(v); changed.add(v); done = True
            for g in occ[v]:
                if g != f and g not in seen: q.append(g); seen.add(g)
            break
        if not done: unf.append(f)
    return changed, unf

def run(base_val, base_fv, base_av, base_eqv, forced, seed, allowed):
    val = list(base_val)
    rnd = random.Random(seed)
    changed, unf = propagate(val, forced, allowed, rnd)
    fv, av, eqv, ok = incr_score(val, changed, base_fv, base_av, base_eqv)
    return ok, val, len(unf), len(changed)

if __name__ == '__main__':
    base = sys.argv[1]
    forced = {}
    for kv in sys.argv[2].split(','):
        k, v = kv.split('='); forced[int(k)] = int(v)
    ntry = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    outp = sys.argv[4] if len(sys.argv) > 4 else None
    val0 = load(base)
    fv0, av0, eqv0 = full_state(val0)
    base_ok = sum(1 for x in eqv0 if x == 0)
    allowed = set(f for f in range(len(facs)) if fv0[f] != 0)
    print("base %d  allowed-nonzero gates %d" % (base_ok, len(allowed)))
    best = (base_ok, None)
    t0 = time.time()
    hist = collections.Counter()
    for s in range(ntry):
        ok, val, nu, nc = run(val0, fv0, av0, eqv0, forced, s, allowed)
        hist[ok] += 1
        if ok > best[0]:
            best = (ok, val)
            print("  seed %d -> %d (unfixable %d, changed %d)" % (s, ok, nu, nc), flush=True)
    print("best %d over %d tries  %.1fs" % (best[0], ntry, time.time()-t0))
    print("score histogram:", sorted(hist.items(), reverse=True)[:10])
    if outp and best[1]:
        json.dump({('x_%d'%i): best[1][i] for i in range(NV)}, open(outp, 'w'))
        print("wrote", outp)
