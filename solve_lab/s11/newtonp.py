"""NEWTON MOD P WITH THE TRUE JACOBIAN.

Every earlier system filtered variables by  f(u+2)-f(u) == 2*(f(u+1)-f(u))  -- i.e. it kept
only variables entering LINEARLY.  That silently discards every quadratic variable, which is
exactly where the mod-p leverage lives (perm.py: 6-30 such variables per failing equation).

Here we use the exact symbolic partial derivative instead, so quadratic variables are included,
and run Newton over GF(p) on the whole region:

        J(x) . delta  =  -residual   (mod p),      x <- x + delta,   iterate.

For a polynomial system this is Hensel/Newton and converges when J is surjective mod p.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip14 import gf_solve
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def datom(a, t, v):
    """exact partial derivative of atom a w.r.t. variable t, at v"""
    s = 0
    for m, c in L.polys[a].items():
        k = m.count(t)
        if not k:
            continue
        term = c * k
        used = False
        for u in m:
            if u == t and not used:
                used = True
                continue
            term *= v[u]
        s += term
    return s


def eq_inner(e, v, AV):
    return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())


def region(v, FAIL, hops=1):
    R = set(FAIL)
    for _ in range(hops):
        vs = set()
        for e in R:
            for a in L.eq_atoms[e][2]:
                vs |= set(L.avars[a])
        for u in vs:
            for a in L.var_atoms[u]:
                R |= set(L.atom2eq.get(a, {}))
    vs = set()
    for e in R:
        for a in L.eq_atoms[e][2]:
            vs |= set(L.avars[a])
    return sorted(R), sorted(vs)


def newton_modp(v, iters=25, hops=1, verbose=True):
    for it in range(iters):
        AV = [atomval(a, v) for a in range(L.NA)]
        FAIL = [e for e in range(L.NEQ) if eq_inner(e, v, AV) != 0]
        badp = [e for e in range(L.NEQ) if eq_inner(e, v, AV) % P != 0]
        if verbose:
            print(f"  it{it}: failing={len(FAIL)}  nonzero mod p={len(badp)}", flush=True)
        if not badp:
            return v, True
        R, VS = region(v, badp, hops)
        # true Jacobian mod p over the region
        idx = {e: i for i, e in enumerate(R)}
        cols = []
        used = []
        for t in VS:
            col = collections.defaultdict(int)
            for a in L.var_atoms[t]:
                d = datom(a, t, v) % P
                if d:
                    for e, c in L.atom2eq.get(a, {}).items():
                        if e in idx:
                            col[idx[e]] = (col[idx[e]] + c * d) % P
            col = {k: x for k, x in col.items() if x}
            if col:
                cols.append(col)
                used.append(t)
        J = [[0] * len(used) for _ in R]
        for j, col in enumerate(cols):
            for i, x in col.items():
                J[i][j] = x
        rhs = [(-eq_inner(e, v, AV)) % P for e in R]
        t0 = time.time()
        d = gf_solve(J, rhs, P)
        if verbose:
            print(f"      region {len(R)} eqs x {len(used)} vars ; GF(p) solve "
                  f"{'OK' if d is not None else 'INCONSISTENT'} ({time.time()-t0:.0f}s)", flush=True)
        if d is None:
            return v, False
        for j, t in enumerate(used):
            v[t] = v[t] + d[j]
    AV = [atomval(a, v) for a in range(L.NA)]
    return v, not [e for e in range(L.NEQ) if eq_inner(e, v, AV) % P != 0]


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    hops = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    v = load_raw(src)
    print(f"=== {os.path.basename(src)}  (true-Jacobian Newton mod p, hops={hops})")
    t0 = time.time()
    v, ok = newton_modp(v, hops=hops)
    AV = [atomval(a, v) for a in range(L.NA)]
    F = [e for e in range(L.NEQ) if eq_inner(e, v, AV) != 0]
    badp = [e for e in range(L.NEQ) if eq_inner(e, v, AV) % P != 0]
    print(f"RESULT: all equations 0 mod p = {ok} ; failing={len(F)} score={L.NEQ-len(F)} ; "
          f"nonzero mod p={len(badp)}  ({time.time()-t0:.0f}s)")
    if ok:
        json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'newtonp_modp.json'), 'w'))
        print("  saved data/newtonp_modp.json  (a state with EVERY equation 0 mod p)")
