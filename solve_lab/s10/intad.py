"""S11 step 76: forward-mode AD over Z, and the handle lift it makes possible.

fidelity.py established that every pivot in the DAG is +-1 (bar one), so forward
evaluation never leaves the integers -- which means forward-mode AD can be run over
Z instead of mod p, with no inverses and no reduction.  That matters because of
§124: a handle enters its check with coefficient d*p, so the mod-p Jacobian reports
ZERO for every handle.  Handles are invisible to all the linear algebra in this lab.
Over Z they are visible, and exact.

Two things follow:
  * `jacZ(u)` gives the exact integer coefficient of every check w.r.t. free input u.
  * if a check is linear in u (all handles are), then c(v + t*e_u) = c(v) + t*gamma
    exactly, so `t = -c/gamma` zeroes it over Z whenever gamma | c.  That is the
    LIFT phase of §124, and it is pure arithmetic -- no search.

`lift(v)` runs it to a fixed point: repeatedly find a nonzero check and a free input
whose exact integer coefficient divides it, verify the move by re-evaluating, and
keep it only if the score does not drop.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import suppfree
P = ad.P


def dpartZ(a, w, v):
    """d(atom a)/d(x_w) over Z.  ad.dpart reduces mod p and therefore reports 0 for
    every p-wire multiplication -- which is exactly where the handles live."""
    s = 0
    for m, c in L.polys[a].items():
        k = m.count(w)
        if k == 0:
            continue
        if k == 1:
            t = c
            for z in m:
                if z != w:
                    t = t * v[z]
            s += t
        else:
            s += 2 * c * v[w]
    return s


def jacZ(u, v, checks):
    """Exact integer d(check)/d(x_u) by one forward pass."""
    dv = {u: 1}
    for t in ad.ORDER:
        a = L.definer[t]
        d = dpartZ(a, t, v)
        if d == 0:
            dv[t] = 0
            continue
        s = 0
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dv.get(w, 0)
            if dw:
                s += dpartZ(a, w, v) * dw
        if s % d:
            dv[t] = None          # not an integer derivative: give up on this path
            continue
        dv[t] = -s // d
    out = {}
    for c in checks:
        s = 0
        ok = True
        for w in L.avars[c]:
            dw = dv.get(w, 0)
            if dw is None:
                ok = False
                break
            if dw:
                s += dpartZ(c, w, v) * dw
        if ok and s:
            out[c] = s
    return out


def lift(v, verbose=True, rounds=6):
    """Absorb every check that an exact integer coefficient can divide."""
    ad.fwd(v, rounds=rounds)
    av = L.all_atom_values(v)
    best = L.NEQ - len(L.failing_eqs(av))
    idx, freelist, vs = suppfree.build(v, modp=False)
    moved = True
    while moved:
        moved = False
        NZ = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
        for c in NZ:
            m = suppfree.atom_supp(c, v, vs, modp=False)
            sup = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]
            for u in sup:
                g = jacZ(u, v, [c]).get(c, 0)
                if not g or av[c] % g:
                    continue
                t = -av[c] // g
                w = list(v)
                w[u] = w[u] + t
                ad.fwd(w, rounds=rounds)
                aw = L.all_atom_values(w)
                s = L.NEQ - len(L.failing_eqs(aw))
                if aw[c] == 0 and s >= best:
                    if verbose:
                        print('    lift: a%d zeroed via x%d (coef %s), score %d -> %d'
                              % (c, u, ('%d' % g)[:24], best, s), flush=True)
                    v, av, best = w, aw, s
                    idx, freelist, vs = suppfree.build(v, modp=False)
                    moved = True
                    break
            if moved:
                break
    return v, best


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else 'mod9118_0.json'
    v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    NZ = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%s: score %d; nonzero checks %s'
          % (src, L.NEQ - len(L.failing_eqs(av)), NZ), flush=True)
    vm = [x % P for x in v]
    idx, freelist, vs = suppfree.build(v, modp=False)
    from fwdad import jac_column
    for c in NZ:
        m = suppfree.atom_supp(c, v, vs, modp=False)
        sup = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]
        modp = jac_column_modp = None
        hidden = []
        for u in sup:
            gp = jac_column(u, v, vm, [c]).get(c, 0)
            gz = jacZ(u, v, [c]).get(c, 0)
            if gp == 0 and gz:
                hidden.append((u, gz))
        print('  a%-6d support %d;  HANDLES invisible mod p: %d'
              % (c, len(sup), len(hidden)), flush=True)
        for u, gz in hidden[:8]:
            q, r = divmod(gz, P)
            print('       x%-6d exact coef = %s*p%s   divides the residual? %s'
                  % (u, q, '' if r == 0 else ' + %d' % r, av[c] % gz == 0))
    v2, s2 = lift(v)
    print('after the lift: score %d' % s2)
    if s2 > L.NEQ - len(L.failing_eqs(av)):
        T.save(v2, os.path.join(HERE, 'lift_%d.json' % s2))
        print('saved lift_%d.json' % s2)
