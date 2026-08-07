"""Exact integer linear system solver via column Hermite normal form."""


def col_hnf(A):
    """A: list of rows (m x n).  Returns (H, U) with A*U = H, U unimodular."""
    m = len(A); n = len(A[0]) if m else 0
    H = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def colop(c1, c2, a, b, c, d):
        """[col c1, col c2] <- [a*c1 + b*c2, c*c1 + d*c2]  (det ad-bc = +-1)"""
        for r in range(m):
            x, y = H[r][c1], H[r][c2]
            H[r][c1] = a * x + b * y
            H[r][c2] = c * x + d * y
        for r in range(n):
            x, y = U[r][c1], U[r][c2]
            U[r][c1] = a * x + b * y
            U[r][c2] = c * x + d * y

    piv = 0
    pivots = []
    for r in range(m):
        if piv >= n:
            break
        # find a nonzero entry in row r among columns >= piv
        nz = [c for c in range(piv, n) if H[r][c] != 0]
        if not nz:
            pivots.append(None)
            continue
        # gcd-reduce columns piv..n-1 in row r into column piv
        for c in nz:
            if c == piv:
                continue
            while H[r][c] != 0:
                if H[r][piv] == 0:
                    colop(piv, c, 0, 1, 1, 0)
                    continue
                q = H[r][c] // H[r][piv]
                # col c -= q * col piv
                colop(piv, c, 1, -q, 0, 1)
                if H[r][c] != 0:
                    colop(piv, c, 0, 1, 1, 0)
        pivots.append(piv)
        piv += 1
    return H, U, pivots


def solve_int(A, b):
    """Return an integer x with A x = b, or None."""
    m = len(A)
    if m == 0:
        return []
    n = len(A[0])
    H, U, pivots = col_hnf(A)
    y = [0] * n
    rhs = b[:]
    for r in range(m):
        p = pivots[r] if r < len(pivots) else None
        s = sum(H[r][c] * y[c] for c in range(n))
        need = rhs[r] - s
        if p is None:
            if need != 0:
                return None
            continue
        if H[r][p] == 0:
            if need != 0:
                return None
            continue
        if need % H[r][p] != 0:
            return None
        y[p] = need // H[r][p]
    # verify
    x = [sum(U[i][c] * y[c] for c in range(n)) for i in range(n)]
    for r in range(m):
        if sum(A[r][c] * x[c] for c in range(n)) != b[r]:
            return None
    return x


if __name__ == '__main__':
    A = [[2, 4], [3, 9]]
    print(solve_int(A, [6, 12]))
    print(solve_int([[2, 4]], [3]))
