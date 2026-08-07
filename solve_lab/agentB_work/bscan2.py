"""Agent B: for every free input, measure the derivative of the residual congruences
and the collateral (newly-violated gates).  Linear-algebra-ready."""
import os, sys, json, pickle, collections, time
os.environ.setdefault('ORIENT', 'orient7.pkl')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentB_work')
import beval as E

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
Q = 7376877 * P
NV = 38748
facs, atoms, eqs = E.facs, E.atoms, E.eqs
occ = collections.defaultdict(list)
for i, p in enumerate(facs):
    vs = set()
    for m in p: vs.update(m)
    for v in vs: occ[v].append(i)

v0 = E.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fv0 = [E.fval(f, v0) for f in range(len(facs))]
base_nonzero = set(f for f in range(len(facs)) if fv0[f] != 0)
freeval0 = {v: v0[v] for v in E.free}

def targets(val):
    return (val[28730] % P, val[8731] % P, val[9118] % P,
            (val[7068]-val[2099]) % Q, val[7927] % P,
            (val[15324]-val[37254]) % 8481759, val[579] % P)

T0 = targets(v0)

def probe(v, d):
    fvl = dict(freeval0); fvl[v] = fvl[v] + d
    val, nd, nf2 = E.forward(fvl, default=v0)
    changed = [u for u in range(NV) if val[u] != v0[u]]
    dirty = set()
    for u in changed: dirty.update(occ[u])
    newbad = [f for f in dirty if f not in base_nonzero and E.fval(f, val) != 0]
    return targets(val), newbad, len(changed), nd

def main():
    out = []
    t0 = time.time()
    for n, v in enumerate(E.free):
        T1, nb1, nc1, nd1 = probe(v, 1)
        d1 = tuple((T1[i] - T0[i]) for i in range(7))
        if any(d1) or nb1:
            T2, nb2, nc2, nd2 = probe(v, 2)
            d2 = tuple((T2[i] - T0[i]) for i in range(7))
            lin = all((2*d1[i] - d2[i]) % m == 0 for i, m in
                      enumerate([P, P, P, Q, P, 8481759, P]))
            out.append((v, d1, sorted(nb1), nc1, nd1, lin))
        if (n+1) % 500 == 0:
            print("  %d/%d  hits=%d  %.1fs" % (n+1, len(E.free), len(out), time.time()-t0), flush=True)
    print("movers:", len(out))
    pickle.dump({'T0': T0, 'rows': out}, open('/home/user/integer_solver/solve_lab/agentB_work/scan2.pkl','wb'), -1)
    names = ['A', 'B', 'E', 'D', 'G1', 'G2', 'G3']
    for v, d1, nb, nc, nd, lin in out:
        act = [names[i] for i in range(7) if d1[i]]
        if act:
            print("  x%-6d moves %-22s collateral=%d lin=%s nchanged=%d" % (v, ','.join(act), len(nb), lin, nc))

if __name__ == '__main__':
    main()
