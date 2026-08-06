"""jm step 8: minimum-weight collateral over F_p, and an EXPLICIT direction.

  * rank of the check matrix, and confirm target in its rowspan
  * search for the minimum-weight set D of checks that must break
    (T not in rowspan(A \\ D)); try |D| = 1, then 2, then the greedy set
  * build a concrete c in ker(A \\ D) with Tc = 1 and dump it
"""
import os, sys, json, time, collections, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L
from jm_07_lin import load, weight, reduce_by, add_to, FREEATOMS
P = J.P
OUT = '/home/user/integer_solver/solve_lab/s10/jm_dirs.jsonl'


def build():
    cols, dR, dA1 = load()
    U = sorted(cols)
    ui = {u: i for i, u in enumerate(U)}
    rows = collections.defaultdict(dict)
    for u, c in cols.items():
        for a, d in c.items():
            rows[a][u] = d
    rows = {a: r for a, r in rows.items() if a not in FREEATOMS}
    return U, ui, rows, cols, dR, dA1


def span_basis(rowlist):
    basis = []
    for r in rowlist:
        nb = add_to(basis, r)
        if nb is not None:
            basis = nb
    return basis


def in_span(T, basis):
    return not reduce_by(T, basis)


def kernel_vec(rowlist, T, U):
    """find c with A c = 0 for rows in rowlist and T c = 1, by Gaussian
    elimination over the columns.  Returns dict u->val or None."""
    n = len(U)
    ui = {u: i for i, u in enumerate(U)}
    # dense-ish sparse gaussian elimination: pivot per row
    piv = {}                     # col -> reduced row (dict)
    for r in rowlist:
        r = {ui[u]: v % P for u, v in r.items() if v % P}
        for c in sorted(r):
            if c in piv:
                f = r[c]
                for k, v in piv[c].items():
                    nv = (r.get(k, 0) - f * v) % P
                    if nv:
                        r[k] = nv
                    else:
                        r.pop(k, None)
        r = {k: v for k, v in r.items() if v}
        if not r:
            continue
        c0 = min(r)
        inv = pow(r[c0], -1, P)
        piv[c0] = {k: v * inv % P for k, v in r.items()}
    # back-substitute to fully reduced form
    for c in sorted(piv, reverse=True):
        row = piv[c]
        for k in sorted([k for k in row if k != c and k in piv]):
            f = row.get(k, 0)
            if not f:
                continue
            for kk, vv in piv[k].items():
                nv = (row.get(kk, 0) - f * vv) % P
                if nv:
                    row[kk] = nv
                else:
                    row.pop(kk, None)
        piv[c] = row
    freecols = [i for i in range(n) if i not in piv]
    Td = {ui[u]: v % P for u, v in T.items() if v % P}
    # express T on the free columns:  c is determined by its free entries
    #   c[pivot] = -sum_{free k} piv[pivot][k]*c[k]
    coef = collections.defaultdict(int)
    for i, val in Td.items():
        if i in piv:
            for k, v in piv[i].items():
                if k != i:
                    coef[k] = (coef[k] - val * v) % P
        else:
            coef[i] = (coef[i] + val) % P
    coef = {k: v for k, v in coef.items() if v % P and k in set(freecols)}
    if not coef:
        return None
    k0 = min(coef)
    c = {k0: pow(coef[k0], -1, P)}
    full = dict(c)
    for pc, row in piv.items():
        s = 0
        for k, v in row.items():
            if k != pc and k in c:
                s += v * c[k]
        if s % P:
            full[pc] = (-s) % P
    return {U[i]: v % P for i, v in full.items() if v % P}


if __name__ == '__main__':
    t0 = time.time()
    U, ui, rows, cols, dR, dA1 = build()
    W = {a: weight(a) for a in rows}
    allrows = sorted(rows)
    print(f'{len(allrows)} rows, {len(U)} cols', flush=True)
    B = span_basis([rows[a] for a in allrows])
    print(f'rank(A) = {len(B)}   kernel dim = {len(U)-len(B)}  '
          f'({time.time()-t0:.0f}s)', flush=True)

    TARGETS = {'C1': {u: dR[u][0] for u in U if dR[u][0]},
               'C2': {u: dA1[u] for u in U if dA1[u]},
               'BOTH_sum': None}
    res = {}
    f = open(OUT, 'a')
    for nm in ('C1', 'C2'):
        T = TARGETS[nm]
        print(f'\n### {nm}: T in rowspan(A)? {in_span(T, B)}', flush=True)
        # --- minimum weight D ---
        best = None
        cand = sorted(allrows, key=lambda a: (W[a], a))
        for a in cand:
            if W[a] > 6:
                break
            Ba = span_basis([rows[b] for b in allrows if b != a])
            if not in_span(T, Ba):
                best = ([a], W[a])
                print(f'  |D|=1 works: break only a{a} (weight {W[a]})')
                break
        if best is None:
            light = [a for a in cand if W[a] <= 3][:70]
            print(f'  no single row; trying pairs among {len(light)} light rows',
                  flush=True)
            for a, b in itertools.combinations(light, 2):
                Ba = span_basis([rows[x] for x in allrows if x not in (a, b)])
                if not in_span(T, Ba):
                    best = ([a, b], W[a] + W[b])
                    print(f'  |D|=2 works: break a{a},a{b} '
                          f'(weight {W[a]+W[b]})')
                    break
        if best is None:
            best = (json.load(open(f'/home/user/integer_solver/solve_lab/s10/'
                                   f'jm_lin_{nm}.json'))['forced'], None)
            best = (best[0], sum(W[a] for a in best[0]))
            print(f'  falling back to greedy D = {best[0]} weight {best[1]}')
        D = set(best[0])
        keep = [rows[a] for a in allrows if a not in D]
        c = kernel_vec(keep, T, U)
        print(f'  direction: {"FOUND" if c else "NONE"}'
              + (f' support {len(c)} inputs' if c else ''), flush=True)
        if c:
            nzrows = [a for a in allrows
                      if sum(rows[a].get(u, 0) * v for u, v in c.items()) % P]
            eqs = set()
            for a in nzrows:
                eqs |= set(L.atom2eq.get(a, {})) - J.E12
            print(f'  predicted broken checks {nzrows} -> {len(eqs)} equations')
            f.write(json.dumps({'target': nm, 'D': sorted(D),
                                'c': {str(k): str(v) for k, v in c.items()},
                                'pred_rows': nzrows,
                                'pred_eqs': sorted(eqs)}) + '\n')
            f.flush()
        res[nm] = best
    f.close()
    print(f'\ndone ({time.time()-t0:.0f}s)')
