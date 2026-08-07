"""Agent B: forward evaluation of the recovered circuit (model5).

Try to drive every gate factor to zero by integer unit propagation.
Deterministic gates first; boolean/free choices are deferred and then guessed.
"""
import pickle, collections, sys, time, json

W = '/home/user/integer_solver/solve_lab/agentB_work/'

def load():
    M = pickle.load(open(W+'model5.pkl','rb'))
    return M['facs'], M['atoms'], M['eqs']

facs, atoms, eqs = load()
NV = 38748

fvars = []
for p in facs:
    s = set()
    for m in p: s.update(m)
    fvars.append(s)
occ = collections.defaultdict(list)
for i, s in enumerate(fvars):
    for v in s: occ[v].append(i)

def reduce_fac(p, val):
    """-> (const, {var:coef} linear part, set of vars in nonlinear monomials)"""
    k = 0; lin = collections.defaultdict(int); nl = set()
    for m, c in p.items():
        cc = c; rem = []
        for v in m:
            x = val.get(v)
            if x is None: rem.append(v)
            else:
                cc *= x
                if cc == 0: break
        if cc == 0: continue
        if not rem: k += cc
        elif len(rem) == 1: lin[rem[0]] += cc
        else: nl.update(rem)
    return k, {v: c for v, c in lin.items() if c}, nl

def score_assignment(val):
    fv = [None]*len(facs)
    for i, p in enumerate(facs):
        s = 0
        for m, c in p.items():
            t = c
            for v in m: t *= val.get(v, 0)
            s += t
        fv[i] = s
    nz = sum(1 for x in fv if x)
    ok = 0; failing = []
    for e, (sc, L, k) in enumerate(eqs):
        s = 0
        for c, a in L:
            t = c
            for f in atoms[a]:
                t *= fv[f]
                if t == 0: break
            s += t
        if s == 0: ok += 1
        else: failing.append(e)
    return ok, failing, nz, fv

def propagate(val, allow=None, verbose=True):
    """unit-propagate all factors (target: every factor == 0)."""
    q = collections.deque(range(len(facs)))
    inq = set(q)
    t0 = time.time(); n = 0
    while q:
        i = q.popleft(); inq.discard(i); n += 1
        k, lin, nl = reduce_fac(facs[i], val)
        if nl or len(lin) != 1: continue
        v, c = next(iter(lin.items()))
        if (-k) % c: continue
        x = (-k)//c
        if v in val: continue
        val[v] = x
        for j in occ[v]:
            if j not in inq: q.append(j); inq.add(j)
    if verbose:
        print("  propagate: %d steps, assigned=%d, %.1fs" % (n, len(val), time.time()-t0), flush=True)
    return val

def main():
    val = {}
    propagate(val)
    ok, failing, nz, fv = score_assignment(val)
    print("after pure propagation: assigned=%d  nonzero gates=%d  eqs ok=%d/%d" % (len(val), nz, ok, len(eqs)))
    # remaining unknown vars
    unk = [v for v in range(NV) if v not in val]
    print("unknown vars:", len(unk))
    # what blocks?  histogram of unresolved factors
    blocked = collections.Counter()
    for i, p in enumerate(facs):
        k, lin, nl = reduce_fac(p, val)
        if nl: blocked['nonlinear:%d' % min(len(nl),4)] += 1
        elif len(lin) == 0:
            blocked['closed_ok' if k == 0 else 'closed_VIOLATED'] += 1
        elif len(lin) == 1:
            v, c = next(iter(lin.items()))
            blocked['unit_nondivisible' if (-k) % c else 'unit_open'] += 1
        else:
            blocked['linear:%d' % min(len(lin),4)] += 1
    print("factor status:", dict(blocked))
    json.dump({('x_%d'%i): val.get(i,0) for i in range(NV)}, open(W+'out/fwd0.json','w'))
    pickle.dump(val, open(W+'fwd0.pkl','wb'), -1)
    print("wrote out/fwd0.json")

if __name__ == '__main__':
    main()
