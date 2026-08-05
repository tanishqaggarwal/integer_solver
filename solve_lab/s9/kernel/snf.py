"""Smith Normal Form based exact integer linear solver: solve A z = b over Z."""

def _ident(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def snf(Ain):
    A = [r[:] for r in Ain]
    m = len(A)
    n = len(A[0])
    U = _ident(m)
    V = _ident(n)

    def swap_rows(a, b):
        A[a], A[b] = A[b], A[a]
        U[a], U[b] = U[b], U[a]

    def swap_cols(a, b):
        for r in A:
            r[a], r[b] = r[b], r[a]
        for r in V:
            r[a], r[b] = r[b], r[a]

    def addrow(dst, src, q):
        for j in range(n):
            A[dst][j] -= q * A[src][j]
        for j in range(m):
            U[dst][j] -= q * U[src][j]

    def addcol(dst, src, q):
        for i in range(m):
            A[i][dst] -= q * A[i][src]
        for i in range(n):
            V[i][dst] -= q * V[i][src]

    t = 0
    while t < min(m, n):
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j]:
                    piv = (i, j)
                    break
            if piv:
                break
        if piv is None:
            break
        i, j = piv
        if i != t:
            swap_rows(t, i)
        if j != t:
            swap_cols(t, j)
        while True:
            changed = False
            for i in range(t + 1, m):
                if A[i][t]:
                    q = A[i][t] // A[t][t]
                    addrow(i, t, q)
                    if A[i][t]:
                        swap_rows(t, i)
                        changed = True
            for j in range(t + 1, n):
                if A[t][j]:
                    q = A[t][j] // A[t][t]
                    addcol(j, t, q)
                    if A[t][j]:
                        swap_cols(t, j)
                        changed = True
            if changed:
                continue
            if all(A[i][t] == 0 for i in range(t + 1, m)) and \
               all(A[t][j] == 0 for j in range(t + 1, n)):
                break
        t += 1
    return A, U, V


def solve_int(A, b):
    """Return an integer z with A z = b, or None if no integer solution exists."""
    if not A:
        return []
    m = len(A)
    n = len(A[0])
    D, U, V = snf(A)
    c = [sum(U[i][k] * b[k] for k in range(m)) for i in range(m)]
    y = [0] * n
    for i in range(min(m, n)):
        d = D[i][i]
        if d == 0:
            if c[i] != 0:
                return None
        else:
            if c[i] % d:
                return None
            y[i] = c[i] // d
    for i in range(min(m, n), m):
        if c[i] != 0:
            return None
    z = [sum(V[r][k] * y[k] for k in range(n)) for r in range(n)]
    for k in range(m):
        if sum(A[k][j] * z[j] for j in range(n)) != b[k]:
            return None
    return z
