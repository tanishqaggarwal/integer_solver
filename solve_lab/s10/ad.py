"""S10 step 33: exact reverse-mode AD of the circuit, mod p.

Every gate atom a defines a variable x_t via a(x) = 0, so
    dx_t = -sum_{w != t} (da/dx_w) / (da/dx_t) * dx_w.
Propagating adjoints in reverse topological order gives d(check)/d(free input)
mod p for every free input, at a cost of one pass per check.

Crucially, mod p every monomial containing a p-wire variable dies, so the
circuit is far simpler mod p than it is over Z.  This module also measures how
nonlinear the check map actually is mod p.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
topo = list(L.topo)
_tset = set(topo)
cyc = [x for x in definer if x not in _tset]
ORDER = cyc + topo               # forward order used by fwd()
FREE = [u for u in range(L.NVARS) if u not in definer]
FREESET = set(FREE)


def dpart(a, w, vm):
    """d(atom a)/d(x_w) mod p, evaluated at vm (list of residues)."""
    s = 0
    for m, c in L.polys[a].items():
        k = m.count(w)
        if k == 0:
            continue
        if k == 1:
            t = c
            for z in m:
                if z != w:
                    t = t * vm[z] % P
            s += t
        else:                                  # w^2
            s += 2 * c * vm[w]
    return s % P


def grad(check, vm):
    """dict free_input -> d(check)/d(free input) mod p."""
    lam = collections.defaultdict(int)
    for w in L.avars[check]:
        lam[w] = (lam[w] + dpart(check, w, vm)) % P
    for t in reversed(ORDER):
        lt = lam.get(t, 0)
        if not lt:
            continue
        a = definer[t]
        d = dpart(a, t, vm)
        if d % P == 0:
            continue
        inv = pow(d, -1, P)
        f = -lt * inv % P
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dpart(a, w, vm)
            if dw:
                lam[w] = (lam[w] + f * dw) % P
        lam[t] = 0
    return {u: lam[u] % P for u in FREESET if lam.get(u, 0) % P}


def fwd(v, rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None:
                v[u] = nv
    return v


if __name__ == '__main__':
    v = L.load(os.path.join(HERE, 'forward_state.json'))
    vm = [x % P for x in v]
    av = L.all_atom_values(v)
    CHECKS = [a for a in range(L.NA) if av[a] and a not in atom_out]
    print('failing checks:', CHECKS)

    for c in CHECKS:
        g = grad(c, vm)
        print(f'\ncheck a{c}: residual mod p = {av[c] % P}')
        print(f'   gradient support: {len(g)} free inputs (of {len(FREE)})')
        top = sorted(g.items(), key=lambda kv: -abs(kv[1]))[:6]
        print(f'   sample: {[(f"x_{u}", str(d)[:16]) for u, d in list(g.items())[:6]]}')

        # ---- validate against an exact finite difference -------------------
        ok = bad = 0
        for u in list(g)[:5]:
            w = list(v)
            w[u] = w[u] + 1
            fwd(w)
            got = (L.evalpoly(L.polys[c], w) - av[c]) % P if hasattr(L, 'evalpoly') else None
            if got is None:
                break
            if got == g[u] % P:
                ok += 1
            else:
                bad += 1
                if bad <= 2:
                    print(f'      MISMATCH x_{u}: AD={g[u]} FD={got}')
        print(f'   finite-difference validation: {ok} match, {bad} mismatch')
