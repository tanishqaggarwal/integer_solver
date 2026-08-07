#!/usr/bin/env python3
"""Option 3 measured: CP encoding over branch choices (OR-Tools CP-SAT), same shape.
CP-SAT domains are int64, so this can only be posed directly for primes up to ~31 bits;
beyond that the field elements do not fit and limb decomposition (= bit-blasting) is forced."""
import sys, time, json
from ortools.sat.python import cp_model
import sibling, witness

def run(m, tmo=300, pin=False, workers=1):
    d = sibling.instance(m); witness.witness(d, d['k'])
    p, n, lad, T, G = d['p'], d['n'], d['lad'], d['T'], d['G']
    assert p * p < 2 ** 62
    mo = cp_model.CpModel()
    def V(nm, hi=None): return mo.NewIntVar(0, (hi if hi is not None else p - 1), nm)
    def mulmod(a, b, nm):
        pr = mo.NewIntVar(0, p * p, nm + '_p'); mo.AddMultiplicationEquality(pr, [a, b])
        r = V(nm); q = V(nm + '_q'); mo.Add(pr == p * q + r)
        return r
    def red(terms, const, nm, lo=4):
        e = sum(terms) + const + lo * p
        r = V(nm); q = mo.NewIntVar(0, 2 * lo, nm + '_r'); mo.Add(e == p * q + r)
        return r
    ax = mo.NewConstant(G[0]); ay = mo.NewConstant(G[1])
    bits = []
    for i in range(1, n):
        bx, by = lad[i]
        bi = mo.NewBoolVar('b%d' % i); bits.append(bi)
        dd = red([-ax], bx, 'd%d' % i); nn = red([-ay], by, 'n%d' % i)
        mo.Add(dd != 0)
        d2 = mulmod(dd, dd, 'd2_%d' % i); n2 = mulmod(nn, nn, 'n2_%d' % i)
        sx = V('sx%d' % i); sy = V('sy%d' % i)
        t = red([sx, ax], bx, 't%d' % i)
        mo.Add(mulmod(t, d2, 'td%d' % i) == n2)
        u = red([sy, ay], 0, 'u%d' % i); vv = red([ax, -sx], 0, 'v%d' % i)
        mo.Add(mulmod(u, dd, 'ud%d' % i) == mulmod(nn, vv, 'nv%d' % i))
        nax = V('ax%d' % i); nay = V('ay%d' % i)
        mo.Add(nax == sx).OnlyEnforceIf(bi); mo.Add(nax == ax).OnlyEnforceIf(bi.Not())
        mo.Add(nay == sy).OnlyEnforceIf(bi); mo.Add(nay == ay).OnlyEnforceIf(bi.Not())
        ax, ay = nax, nay
    mo.Add(ax == T[0]); mo.Add(ay == T[1])
    if pin:
        for i, bi in enumerate(bits, 1): mo.Add(bi == ((d['k'] >> i) & 1))
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = tmo
    sv.parameters.num_search_workers = workers
    t0 = time.time()
    st = sv.Solve(mo)
    return dict(m=m, n=n, p=p, pin=pin, status=sv.StatusName(st),
                time=round(time.time() - t0, 2),
                conflicts=sv.NumConflicts(), branches=sv.NumBranches(),
                booleans=sv.ResponseProto().solution_info if False else None,
                wall=round(sv.WallTime(), 2))

if __name__ == '__main__':
    out = {}
    for m in (8, 10, 12, 14, 16, 20, 24, 28, 31):
        for pin in (True, False):
            r = run(m, tmo=int(sys.argv[1]) if len(sys.argv) > 1 else 300, pin=pin)
            out['m%d_%s' % (m, 'pin' if pin else 'free')] = r
            json.dump(out, open('runs/cpsat.json', 'w'), indent=1)
            print('m=%2d pin=%-5s %-12s %7.1fs  conflicts=%d branches=%d'
                  % (m, pin, r['status'], r['time'], r['conflicts'], r['branches']), flush=True)
