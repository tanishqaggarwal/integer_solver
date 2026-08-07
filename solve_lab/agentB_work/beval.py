"""Agent B: forward evaluation of the oriented circuit from free-input values."""
import pickle, collections, json, sys, time

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
import os
O = pickle.load(open(W+os.environ.get('ORIENT','orient2.pkl'),'rb'))
matchF, order, free, assertions = O['matchF'], O['order'], O['free'], O['assertions']
NV = 38748

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
    if a == 0: return None if k else 0
    if (-k) % a: return None
    return (-k)//a

def forward(freeval, default=None):
    """freeval: dict var->int for the free inputs.  `default`: fallback values for
    variables whose defining gate degenerates (coefficient 0)."""
    val = [0]*NV if default is None else list(default)
    for v in range(NV):
        if default is None: break
    for v, x in freeval.items(): val[v] = x
    nd = 0; nfree2 = 0
    for f in order:
        v = matchF[f]
        x = solve_for(f, v, val)
        if x is None:
            nd += 1
            x = default[v] if default is not None else 0
        elif x == 0 and _degenerate(f, v, val):
            nfree2 += 1
            x = default[v] if default is not None else 0
        val[v] = x
    return val, nd, nfree2

def _degenerate(f, v, val):
    a = 0
    for m, c in facs[f].items():
        if v in m:
            z = c
            for u in m:
                if u == v: continue
                z *= val[u]
                if z == 0: break
            a += z
    return a == 0

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

if __name__ == '__main__':
    base = sys.argv[1]
    v0 = load(base)
    ok0, fail0, fv0 = score(v0)
    print("base score", ok0)
    freeval = {v: v0[v] for v in free}
    t0 = time.time()
    val, nd, nf2 = forward(freeval, default=v0)
    print("forward eval: %d non-integral divisions, %d degenerate gates  %.1fs" % (nd, nf2, time.time()-t0))
    diff = [v for v in range(NV) if val[v] != v0[v]]
    print("vars differing from base:", len(diff))
    ok, fail, fv = score(val)
    print("forward score:", ok, fail[:15])
    nzass = [f for f in assertions if fv[f] != 0]
    print("assertion gates: %d, violated: %d" % (len(assertions), len(nzass)))
    nzass0 = [f for f in assertions if fv0[f] != 0]
    print("  (at base, violated assertions: %d)" % len(nzass0))
    if len(sys.argv) > 2:
        json.dump({('x_%d'%i): val[i] for i in range(NV)}, open(sys.argv[2],'w'))
        print("wrote", sys.argv[2])
