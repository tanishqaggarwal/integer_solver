"""Exact integer linear system solver over Z via column HNF."""

def col_hnf(M, n):
    """M: list of rows (lists of ints), n columns.  Column-reduce with unimodular
    column ops; return (M, U, pivots) with M_orig @ U = M and pivots=[(row,col)]."""
    m = len(M)
    M = [row[:] for row in M]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def colop(c, c0, q):          # col_c -= q * col_c0
        if q == 0:
            return
        for i in range(m):
            M[i][c] -= q * M[i][c0]
        for i in range(n):
            U[i][c] -= q * U[i][c0]

    def swap(c1, c2):
        if c1 == c2:
            return
        for i in range(m):
            M[i][c1], M[i][c2] = M[i][c2], M[i][c1]
        for i in range(n):
            U[i][c1], U[i][c2] = U[i][c2], U[i][c1]

    piv = 0
    pivots = []
    for r in range(m):
        if piv >= n:
            break
        while True:
            nz = [c for c in range(piv, n) if M[r][c] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda c: abs(M[r][c]))
            c0 = nz[0]
            for c in nz[1:]:
                colop(c, c0, M[r][c] // M[r][c0])
        nz = [c for c in range(piv, n) if M[r][c] != 0]
        if nz:
            swap(piv, nz[0])
            pivots.append((r, piv))
            piv += 1
    return M, U, pivots


def solve_int(Mo, b, n):
    """Find integer x with Mo x = b.  Returns (x, nullbasis) or (None, reason)."""
    M, U, pivots = col_hnf(Mo, n)
    m = len(M)
    y = [0] * n
    pr = {r: c for r, c in pivots}
    for r in range(m):
        s = sum(M[r][c] * y[c] for c in range(n) if y[c])
        need = b[r] - s
        if r in pr:
            c = pr[r]
            if need % M[r][c] != 0:
                return None, f'row {r}: {need} not divisible by {M[r][c]}'
            y[c] = need // M[r][c]
        else:
            if need != 0:
                return None, f'row {r}: inconsistent ({need} != 0)'
    x = [sum(U[i][c] * y[c] for c in range(n)) for i in range(n)]
    rank = len(pivots)
    null = [[U[i][c] for i in range(n)] for c in range(rank, n)]
    return x, null
