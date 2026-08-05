"""Fast max-fixable-subset over F_P: eliminate the must-rows once, then rank-test
all subsets of the (few) candidate rows."""
import sys, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

P = lib.P


def rref_modp(rows, n):
    """rows: list of length-(n+1) lists (augmented). In place RREF mod P.
    Returns (rows, pivots)."""
    rows = [[c % P for c in r] for r in rows]
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pr = i; break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        inv = pow(rows[r][c], P - 2, P)
        rows[r] = [(x * inv) % P for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(a - f * b) % P for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows, piv


def max_subset(must_rows, cand_rows, n):
    """must_rows: list of length-(n+1) augmented rows that must hold (rhs included).
    cand_rows: list of (eqid, augmented row).  Returns (k, best_subset_ids)."""
    M, piv = rref_modp(must_rows, n) if must_rows else ([], [])
    # reduce candidate rows against must
    red = []
    for eid, row in cand_rows:
        r = [c % P for c in row]
        for i, c in enumerate(piv):
            if r[c]:
                f = r[c]
                r = [(a - f * b) % P for a, b in zip(r, M[i])]
        red.append((eid, r))
    m = len(red)
    best = (0, ())
    for k in range(m, 0, -1):
        for sub in itertools.combinations(range(m), k):
            rows = [red[i][1] for i in sub]
            R, pv = rref_modp(rows, n)
            ok = True
            for r in R:
                if all(r[c] == 0 for c in range(n)) and r[n] != 0:
                    ok = False; break
            if ok:
                return k, tuple(red[i][0] for i in sub)
        # continue to smaller k
    return best
