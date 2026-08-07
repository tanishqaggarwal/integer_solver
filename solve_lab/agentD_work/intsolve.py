"""Integer linear system solver via column Hermite normal form."""


def colhnf(A):
    """Return (H, V) with A*V = H in column echelon form, V unimodular."""
    A = [row[:] for row in A]
    m = len(A)
    n = len(A[0]) if m else 0
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def addcol(j, k, c):
        for r in A:
            r[j] += c * r[k]
        for r in V:
            r[j] += c * r[k]

    def swapcol(j, k):
        for r in A:
            r[j], r[k] = r[k], r[j]
        for r in V:
            r[j], r[k] = r[k], r[j]

    piv = 0
    for i in range(m):
        if piv >= n:
            break
        while True:
            nz = [j for j in range(piv, n) if A[i][j] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda j: abs(A[i][j]))
            j0 = nz[0]
            for j in nz[1:]:
                q = A[i][j] // A[i][j0]
                addcol(j, j0, -q)
        nz = [j for j in range(piv, n) if A[i][j] != 0]
        if not nz:
            continue
        if nz[0] != piv:
            swapcol(piv, nz[0])
        piv += 1
    return A, V


def solve_int(A, b):
    """Solve A x = b over Z; return x (list) or None."""
    m = len(A)
    n = len(A[0]) if m else 0
    if n == 0:
        return None if any(b) else []
    H, V = colhnf(A)
    y = [0] * n
    r = list(b)
    col = 0
    for i in range(m):
        j = None
        for jj in range(col, n):
            if H[i][jj] != 0 and all(H[k][jj] == 0 for k in range(i)):
                j = jj
                break
        if j is None:
            continue
        if r[i] % H[i][j]:
            return None
        c = r[i] // H[i][j]
        y[j] = c
        for k in range(m):
            r[k] -= c * H[k][j]
        col = j + 1
    if any(x != 0 for x in r):
        return None
    x = [sum(V[a][c] * y[c] for c in range(n)) for a in range(n)]
    for i in range(m):
        if sum(A[i][k] * x[k] for k in range(n)) != b[i]:
            return None
    return x
