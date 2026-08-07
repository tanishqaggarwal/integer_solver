"""Exact integer linear system solver: A x = b over Z.

Column Hermite normal form with a modular pre-filter.  The pre-filter is what
makes the support search affordable: solvability over Z implies solvability over
F_q for every prime q, and the reduced systems are tiny.
"""
from math import gcd

SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 1000003)
P256 = 2**256 - 2**32 - 977
# prime factors of 7376877 (a modulus that appears in the residual rows) plus p
FILTER_PRIMES = SMALL_PRIMES + (819653, P256)


def _rank_mod(rows, q):
    """Rank of the matrix over Z/q for PRIME q (rows already reduced)."""
    m = len(rows)
    if m == 0:
        return 0
    n = len(rows[0])
    A = [[x % q for x in r] for r in rows]
    rank = 0
    for c in range(n):
        piv = None
        for r in range(rank, m):
            if A[r][c] % q:
                piv = r
                break
        if piv is None:
            continue
        A[rank], A[piv] = A[piv], A[rank]
        inv = pow(A[rank][c], -1, q)
        A[rank] = [x * inv % q for x in A[rank]]
        for r in range(m):
            if r != rank and A[r][c]:
                f = A[r][c]
                A[r] = [(A[r][k] - f * A[rank][k]) % q for k in range(n)]
        rank += 1
        if rank == m:
            break
    return rank


def _mod_consistent(A, b, q):
    """Necessary condition: A x = b solvable over F_q."""
    ra = _rank_mod(A, q)
    rab = _rank_mod([A[i] + [b[i]] for i in range(len(A))], q)
    return ra == rab


def col_hnf(A):
    """A: m x n list of rows.  Returns (H, U, pivots) with A*U = H."""
    m = len(A)
    n = len(A[0]) if m else 0
    H = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def addcol(dst, src, q):
        """col_dst -= q * col_src"""
        if q == 0:
            return
        for r in range(m):
            H[r][dst] -= q * H[r][src]
        for r in range(n):
            U[r][dst] -= q * U[r][src]

    def swapcol(c1, c2):
        for r in range(m):
            H[r][c1], H[r][c2] = H[r][c2], H[r][c1]
        for r in range(n):
            U[r][c1], U[r][c2] = U[r][c2], U[r][c1]

    piv = 0
    pivots = []
    for r in range(m):
        if piv >= n:
            pivots.append(None)
            continue
        while True:
            nz = [c for c in range(piv, n) if H[r][c] != 0]
            if not nz:
                break
            # move the smallest-magnitude nonzero into the pivot column
            c0 = min(nz, key=lambda c: abs(H[r][c]))
            if c0 != piv:
                swapcol(piv, c0)
            done = True
            for c in range(piv + 1, n):
                if H[r][c] != 0:
                    addcol(c, piv, H[r][c] // H[r][piv])
                    if H[r][c] != 0:
                        done = False
            if done:
                break
        if H[r][piv] != 0:
            pivots.append(piv)
            piv += 1
        else:
            pivots.append(None)
    return H, U, pivots


def solve_int(A, b, use_filter=True):
    """Return an integer x with A x = b, or None."""
    m = len(A)
    if m == 0:
        return []
    n = len(A[0])
    if use_filter:
        for q in FILTER_PRIMES:
            if not _mod_consistent(A, b, q):
                return None
    H, U, pivots = col_hnf(A)
    y = [0] * n
    for r in range(m):
        s = sum(H[r][c] * y[c] for c in range(n))
        need = b[r] - s
        p = pivots[r]
        if p is None:
            if need != 0:
                return None
            continue
        if need % H[r][p] != 0:
            return None
        y[p] = need // H[r][p]
    x = [sum(U[i][c] * y[c] for c in range(n)) for i in range(n)]
    for r in range(m):
        if sum(A[r][c] * x[c] for c in range(n)) != b[r]:
            return None
    return x


if __name__ == '__main__':
    assert solve_int([[2, 4], [3, 9]], [6, 12]) is not None
    assert solve_int([[2, 4]], [3]) is None
    assert solve_int([[2, 3]], [1]) is not None
    x = solve_int([[6, 10, 15]], [1])
    assert x and 6 * x[0] + 10 * x[1] + 15 * x[2] == 1
    print("intsolve self-tests OK")
