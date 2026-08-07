"""Agent B: gate-level repair.  Start from an assignment, force some variables,
then greedily restore every violated gate by adjusting a single variable that
occurs linearly with coefficient +-1, choosing the option that creates the fewest
new violations.  Exact integers throughout.
"""
import pickle, collections, json, sys, time, heapq

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

# linear-solvable candidates per gate: var with a degree-1 monomial coeff +-1, not in any other monomial
cands = []
for p in facs:
    quad = set()
    for m in p:
        if len(m) > 1: quad.update(m)
    cands.append([m[0] for m, c in p.items() if len(m) == 1 and abs(c) == 1 and m[0] not in quad])

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
    """value of v making gate f zero, or None"""
    k = 0; a = 0
    for m, c in facs[f].items():
        if m == (v,):
            a += c; continue
        z = c
        for u in m:
            z *= val[u]
            if z == 0: break
        k += z
    if a == 0: return None
    if (-k) % a: return None
    return (-k)//a

def load(path):
    val = [0]*NV
    for k, v in json.load(open(path)).items(): val[int(k.split('_')[1])] = int(v)
    return val

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

def repair(val, frozen, maxiter=400000, verbose=True):
    bad = set(f for f in range(len(facs)) if fval(f, val) != 0)
    t0 = time.time(); it = 0
    stuck = set()
    while True:
        work = [f for f in bad if f not in stuck]
        if not work: break
        prog = False
        for f in sorted(work, key=lambda f: len(cands[f])):
            it += 1
            if fval(f, val) == 0:
                bad.discard(f); prog = True; continue
            best = None
            for v in cands[f]:
                if v in frozen: continue
                x = solve_for(f, v, val)
                if x is None: continue
                old = val[v]
                if old == x: continue
                val[v] = x
                nb = sum(1 for g in occ[v] if fval(g, val) != 0)
                ob = sum(1 for g in occ[v] if (val.__setitem__(v, old), fval(g, val), val.__setitem__(v, x))[1] != 0)
                val[v] = old
                delta = nb - ob
                if best is None or delta < best[0]:
                    best = (delta, v, x)
            if best is None:
                stuck.add(f); continue
            d, v, x = best
            val[v] = x
            for g in occ[v]:
                if fval(g, val) != 0: bad.add(g)
                else: bad.discard(g)
                stuck.discard(g)
            prog = True
            if it % 500 == 0 and verbose:
                print("   it=%d bad=%d stuck=%d %.1fs" % (it, len(bad), len(stuck), time.time()-t0), flush=True)
            if it > maxiter: return val, bad
        if not prog: break
    return val, bad

if __name__ == '__main__':
    base = sys.argv[1]
    val = load(base)
    ok, fail, fv = score(val)
    print("base score %d, nonzero gates %d" % (ok, sum(1 for x in fv if x)))
