"""jm step 7: over F_p, which check atoms MUST break when a congruence moves?

Rows  = check atoms (excluding the 7 residual + a22231, whose equations are all
        inside the twelve, so they are not collateral)
Cols  = 957 free inputs
Target= dC0 (congruence 1)  or  e_28730 (congruence 2)

For a move direction c the collateral is {i : A_i c != 0}.  Greedily grow the set
S of checks forced to zero, heaviest first, keeping the target reachable
(target not in rowspan(A_S)).  What is left is a linear LOWER BOUND on collateral
-- to be confirmed or refuted by construction.
"""
import os, sys, json, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L
P = J.P
JAC = '/home/user/integer_solver/solve_lab/s10/jm_jac.jsonl'
FREEATOMS = set(J.SEVEN) | {22231}


def load():
    cols, dR, dA1 = {}, {}, {}
    for ln in open(JAC):
        r = json.loads(ln)
        u = r['u']
        cols[u] = {int(k): int(v) for k, v in r['col'].items()}
        dR[u] = tuple(int(x) for x in r['dR'])
        dA1[u] = int(r['dA1'])
    return cols, dR, dA1


def weight(a):
    return len(set(L.atom2eq.get(a, {})) - J.E12)


def rref(rows, n):
    """rows: list of dict col->val.  returns list of (pivot, densedict)"""
    basis = []            # list of (piv, row dict) with row[piv]==1
    for r in rows:
        r = dict(r)
        for piv, br in basis:
            if r.get(piv):
                f = r[piv]
                for k, v in br.items():
                    nv = (r.get(k, 0) - f * v) % P
                    if nv:
                        r[k] = nv
                    else:
                        r.pop(k, None)
        r = {k: v % P for k, v in r.items() if v % P}
        if not r:
            continue
        piv = min(r)
        inv = pow(r[piv], -1, P)
        r = {k: v * inv % P for k, v in r.items()}
        basis.append((piv, r))
    return basis


def reduce_by(r, basis):
    r = dict(r)
    for piv, br in basis:
        if r.get(piv):
            f = r[piv]
            for k, v in br.items():
                nv = (r.get(k, 0) - f * v) % P
                if nv:
                    r[k] = nv
                else:
                    r.pop(k, None)
    return {k: v % P for k, v in r.items() if v % P}


def add_to(basis, r):
    """returns new basis if r independent, else None"""
    rr = reduce_by(r, basis)
    if not rr:
        return None
    piv = min(rr)
    inv = pow(rr[piv], -1, P)
    rr = {k: v * inv % P for k, v in rr.items()}
    return basis + [(piv, rr)]


if __name__ == '__main__':
    t0 = time.time()
    cols, dR, dA1 = load()
    U = sorted(cols)
    print(f'{len(U)} columns loaded ({time.time()-t0:.0f}s)', flush=True)

    # transpose: check -> {u: d}
    rows = collections.defaultdict(dict)
    for u, c in cols.items():
        for a, d in c.items():
            rows[a][u] = d
    rows = {a: r for a, r in rows.items() if a not in FREEATOMS}
    print(f'{len(rows)} check rows with support', flush=True)

    W = {a: weight(a) for a in rows}
    print('heaviest rows:',
          sorted(((w, a) for a, w in W.items()), reverse=True)[:14])
    print('zero-weight rows (all eqs inside the twelve):',
          [a for a in rows if W[a] == 0])

    TARGETS = {'C1 (dC0)': {u: dR[u][0] for u in U if dR[u][0]},
               'C2 (dA1)': {u: dA1[u] for u in U if dA1[u]},
               'RB': {u: dR[u][1] for u in U if dR[u][1]},
               'RC': {u: dR[u][2] for u in U if dR[u][2]}}
    for nm, T in TARGETS.items():
        print(f'\n### target {nm}: support {len(T)} inputs '
              f'{sorted(T)[:12]}{"..." if len(T) > 12 else ""}')
        order = sorted(rows, key=lambda a: (-W[a], a))
        basis, forced, killed = [], [], []
        for a in order:
            nb = add_to(basis, rows[a])
            if nb is None:            # a is already in the span -> zeroing free
                killed.append(a)
                continue
            # would adding row a put T in the span?
            if not reduce_by(T, nb):
                forced.append(a)      # cannot zero this one
            else:
                basis = nb
                killed.append(a)
        cost = sum(W[a] for a in forced)
        eqs = set()
        for a in forced:
            eqs |= set(L.atom2eq.get(a, {})) - J.E12
        print(f'  linear floor: must break {forced} '
              f'(weights {[W[a] for a in forced]}) -> {len(eqs)} distinct eqs')
        print(f'  ({time.time()-t0:.0f}s)', flush=True)
        json.dump({'forced': forced, 'weights': [W[a] for a in forced],
                   'eqs': sorted(eqs)},
                  open(f'/home/user/integer_solver/solve_lab/s10/jm_lin_'
                       f'{nm.split()[0]}.json', 'w'))
