"""S10 step 101: FORWARD-mode AD mod p -> the exact Jacobian of every check
with respect to a set of free inputs, one DAG pass per input.

Reverse mode costs one pass per CHECK (10,792 of them).  Forward mode costs one
pass per FREE INPUT, and we only care about a few hundred.  That makes the whole
check-by-free-input Jacobian affordable, and the repair becomes linear algebra
over F_p instead of a search.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
ORDER = ad.ORDER

def jac_column(u, v, vm, checks):
    """d(check)/d(x_u) mod p for every check, by one forward pass."""
    dv = collections.defaultdict(int)
    dv[u] = 1
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0:
            dv[t] = 0; continue
        s = 0
        for w in L.avars[a]:
            if w == t: continue
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(a, w, vm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    out = {}
    for c in checks:
        s = 0
        for w in L.avars[c]:
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(c, w, vm) * dw
        if s % P: out[c] = s % P
    return out

if __name__ == '__main__':
    v = L.load(os.path.join(HERE, 'mod9118_0.json'))
    vm = [x % P for x in v]
    av = L.all_atom_values(v)
    CHECKS = [a for a in range(L.NA) if a not in atom_out]
    BAD = [a for a in CHECKS if av[a]]
    print(f'{len(CHECKS)} check atoms; failing {BAD}')
    # free inputs that can move either failing check
    U = set()
    for a in BAD:
        U |= set(ad.grad(a, vm))
    U = sorted(U - ad.FREE.__class__(( )) if False else U)
    U = [u for u in U if u not in (2081, 4287)]
    print(f'free inputs able to move a failing check: {len(U)}')
    t0 = time.time()
    cols = {}
    for i, u in enumerate(U):
        cols[u] = jac_column(u, v, vm, CHECKS)
        if i % 25 == 0:
            print(f'  col {i}/{len(U)}  ({time.time()-t0:.0f}s) '
                  f'support {len(cols[u])}', flush=True)
    rows = sorted(set().union(*[set(c) for c in cols.values()]))
    print(f'\nJacobian: {len(rows)} affected checks x {len(U)} free inputs '
          f'({time.time()-t0:.0f}s)')
    print(f'   of the affected checks, currently failing: {[r for r in rows if av[r]]}')
    json.dump({'U': U, 'rows': rows,
               'cols': {str(u): {str(c): str(d) for c, d in cols[u].items()} for u in U},
               'resid': {str(r): str(av[r] % P) for r in rows}},
              open(os.path.join(HERE, 'jac.json'), 'w'))
    print('saved jac.json')
