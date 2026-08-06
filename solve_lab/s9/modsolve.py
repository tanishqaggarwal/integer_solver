"""Sparse Gaussian elimination mod p to solve J.delta = -base."""
import pickle, sys, time, collections
P = 2**256-2**32-977

def solve(J, base, P=P, extra_rows=None):
    rows = collections.defaultdict(dict)      # atom -> {col: val}
    for f, col in J.items():
        for a, dv in col.items():
            x = dv % P
            if x: rows[a][f] = x
    rhs = {a: (-base.get(a, 0)) % P for a in rows}
    if extra_rows:
        for a, (rowdict, r) in extra_rows.items():
            rows[a] = dict(rowdict); rhs[a] = r % P
    order = sorted(rows, key=lambda a: len(rows[a]))
    pivots = {}    # col -> (row, rowdict, rhs)
    t0 = time.time()
    for n, a in enumerate(order):
        row = dict(rows[a]); r = rhs[a]
        # reduce by existing pivots
        while True:
            hit = None
            for c in row:
                if c in pivots: hit = c; break
            if hit is None: break
            pr, prow, prhs = pivots[hit]
            factor = row[hit] * pow(prow[hit], P-2, P) % P
            for c, val in prow.items():
                nv = (row.get(c, 0) - factor*val) % P
                if nv: row[c] = nv
                elif c in row: del row[c]
            r = (r - factor*prhs) % P
        if not row:
            if r: return None, f'INCONSISTENT at atom {a}'
            continue
        c = min(row, key=lambda c: len(str(c)))
        pivots[c] = (a, row, r)
    # back-substitute: assign free (non-pivot) columns = 0
    sol = {}
    for c, (a, row, r) in sorted(pivots.items(), key=lambda kv: -len(kv[1][1])):
        pass
    # solve iteratively: since pivot rows are already reduced against each other, do
    # standard back substitution over the pivot set
    piv_cols = list(pivots)
    changed = True
    val = {c: 0 for c in piv_cols}
    for _ in range(len(piv_cols)+2):
        for c, (a, row, r) in pivots.items():
            s = r
            for cc, vv in row.items():
                if cc == c: continue
                s = (s - vv*val.get(cc, 0)) % P
            val[c] = s * pow(row[c], P-2, P) % P
    # verify
    for a, row in rows.items():
        s = sum(v*val.get(c, 0) for c, v in row.items()) % P
        if s != rhs[a]: return None, f'verify failed at atom {a}'
    return val, f'ok ({len(pivots)} pivots, {time.time()-t0:.1f}s)'

if __name__ == '__main__':
    d = pickle.load(open('jac.pkl','rb')); J = d['J']; base = d['base']
    sol, msg = solve(J, base)
    print(msg)
    if sol is not None:
        nz = {c: v for c, v in sol.items() if v}
        print(f'solution: {len(nz)} nonzero free-input deltas (mod p)')
        for c, v in list(nz.items())[:20]: print(f'   x_{c} += {v}')
        pickle.dump(nz, open('modsol.pkl','wb'))
