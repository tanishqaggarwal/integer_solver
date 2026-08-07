"""int_kernel_columns lifted VERBATIM from kerquad.py (whose module body loads a pkl this
measurement does not need).  Byte-identical function text; no reimplementation."""


def int_kernel_columns(A, n):
    """Saturated integer kernel of the integer matrix A (list of rows, n columns).

    Unimodular COLUMN reduction: keep U = I_n and apply the same column operations to U; the
    columns of U whose image column in A has become zero are a basis of ker_Z(A)."""
    m = len(A)
    cols = [[A[i][j] for i in range(m)] for j in range(n)]          # column-major
    U = [[1 if i == j else 0 for i in range(n)] for j in range(n)]  # U[j] = j-th column of U
    piv = 0
    for r in range(m):
        # find a column >= piv with a nonzero entry in row r
        idx = [j for j in range(piv, n) if cols[j][r] != 0]
        if not idx:
            continue
        # gcd-reduce those columns into one
        j0 = idx[0]
        cols[piv], cols[j0] = cols[j0], cols[piv]
        U[piv], U[j0] = U[j0], U[piv]
        for j in range(piv + 1, n):
            while cols[j][r] != 0:
                a, b = cols[piv][r], cols[j][r]
                if abs(a) > abs(b) or a == 0:
                    cols[piv], cols[j] = cols[j], cols[piv]
                    U[piv], U[j] = U[j], U[piv]
                    a, b = cols[piv][r], cols[j][r]
                q = b // a
                if q:
                    cols[j] = [x - q * y for x, y in zip(cols[j], cols[piv])]
                    U[j] = [x - q * y for x, y in zip(U[j], U[piv])]
                else:
                    break
        piv += 1
        if piv == n:
            break
    ker = [U[j] for j in range(n) if all(x == 0 for x in cols[j])]
    return ker


