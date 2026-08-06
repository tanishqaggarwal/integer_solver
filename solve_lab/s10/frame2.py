"""S10 step 106: rebuild the analysis in the DELIVERED witness's orientation.

The delivered witness breaks 5 gate atoms (22229, 22230, 35758, 35761, 35762).
Detach their output variables -- x_7068, x_28730, x_29854, x_31864, x_642 -- so
they become FREE PARAMETERS and those atoms become CHECKS.  In that frame the
delivered witness is on-manifold, the whole gadget cluster is satisfied, and the
residual is exactly the 7 p-quantisation checks.  Then run the same forward-AD
closure and ask whether THOSE are repairable.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
DETACH = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762}
definer = {t: a for t, a in L.definer.items() if t not in DETACH}
atom_out = {a: o for a, o in L.atom_out.items() if a not in set(DETACH.values())}
ORDER = [t for t in ad.ORDER if t not in DETACH]
FREE = set(t for t in range(L.NVARS) if t not in definer)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
print(f'detached {len(DETACH)} variables; free params {len(FREE)} '
      f'(was {len(ad.FREE)}); checks {len(CHECKS)} (was '
      f'{len([a for a in range(L.NA) if a not in L.atom_out])})')

def fwd(v, rounds=6):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None: v[u] = nv
    return v

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
s_before = score(v)
w = list(v); fwd(w)
print(f'delivered witness score {s_before}; after fwd in THIS frame {score(w)} '
      f'(on-manifold: {score(w) == s_before})')
av = L.all_atom_values(w)
BAD = [a for a in CHECKS if av[a]]
print(f'failing checks in this frame: {BAD}')

def grad(check, vm):
    lam = collections.defaultdict(int)
    for x in L.avars[check]:
        lam[x] = (lam[x] + ad.dpart(check, x, vm)) % P
    for t in reversed(ORDER):
        lt = lam.get(t, 0)
        if not lt: continue
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: continue
        f = -lt * pow(d, -1, P) % P
        for x in L.avars[a]:
            if x == t: continue
            dx = ad.dpart(a, x, vm)
            if dx: lam[x] = (lam[x] + f * dx) % P
        lam[t] = 0
    return {u: lam[u] % P for u in FREE if lam.get(u, 0) % P}

def jac_column(u, vm):
    dv = collections.defaultdict(int); dv[u] = 1
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: dv[t] = 0; continue
        s = 0
        for x in L.avars[a]:
            if x == t: continue
            dx = dv.get(x, 0)
            if dx: s += ad.dpart(a, x, vm) * dx
        dv[t] = (-s % P) * pow(d, -1, P) % P
    out = {}
    for c in CHECKS:
        s = 0
        for x in L.avars[c]:
            dx = dv.get(x, 0)
            if dx: s += ad.dpart(c, x, vm) * dx
        if s % P: out[c] = s % P
    return out

if __name__ == '__main__':
    vm = [x % P for x in w]
    U = set()
    for a in BAD: U |= set(grad(a, vm))
    print(f'free params able to move a failing check: {len(U)}')
    print(f'   detached vars among them: {sorted(U & set(DETACH))}')
    cols, t0 = {}, time.time()
    for it in range(10):
        Us = sorted(U)
        for u in Us:
            if u not in cols: cols[u] = jac_column(u, vm)
        rows = sorted(set().union(*[set(cols[u]) for u in Us]) | set(BAD))
        ri = {c: i for i, c in enumerate(rows)}
        n, m = len(rows), len(Us)
        M = [[0] * (m + 1) for _ in rows]
        for j, u in enumerate(Us):
            for c, d in cols[u].items(): M[ri[c]][j] = d % P
        for c in rows: M[ri[c]][m] = (-av[c]) % P
        rid = list(rows); piv, r_ = [], 0
        for j in range(m):
            k = next((i for i in range(r_, n) if M[i][j]), None)
            if k is None: continue
            M[r_], M[k] = M[k], M[r_]; rid[r_], rid[k] = rid[k], rid[r_]
            inv = pow(M[r_][j], -1, P)
            M[r_] = [x * inv % P for x in M[r_]]
            for i in range(n):
                if i != r_ and M[i][j]:
                    f = M[i][j]
                    M[i] = [(a2 - f * b2) % P for a2, b2 in zip(M[i], M[r_])]
            piv.append(j); r_ += 1
        bad_rows = [rid[i] for i in range(r_, n) if M[i][m]]
        print(f'it{it}: {n} rows x {m} cols  rank {r_}  inconsistent {len(bad_rows)} '
              f'({time.time()-t0:.0f}s)', flush=True)
        if not bad_rows:
            d = [0] * m
            for i, j in enumerate(piv): d[j] = M[i][m]
            sol = {Us[j]: d[j] for j in range(m) if d[j]}
            print(f'  *** CONSISTENT -- solution moves {len(sol)} params')
            json.dump({str(u): str(x) for u, x in sol.items()},
                      open(os.path.join(HERE, 'delta2.json'), 'w'))
            print('  saved delta2.json'); break
        new = set()
        for c in bad_rows[:40]: new |= set(grad(c, vm))
        new -= U
        print(f'  witnesses {bad_rows[:8]}; new cols {len(new)}', flush=True)
        if not new:
            print('  FIXED POINT: inconsistent over every reachable parameter'); break
        U |= new
