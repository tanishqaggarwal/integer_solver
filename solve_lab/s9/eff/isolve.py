"""Fast exact integer linear solver:  A z = b over Z, via column Hermite form.

Much faster than the full-SNF solve_int in kernel/snf.py when the entries are large
(the SNF there alternates row/column reduction and suffers coefficient explosion).
Cross-checked against snf.solve_int in selftest().
"""


def solve_int2(A, b):
    """Return an integer z with A z = b, or None if there is no integer solution."""
    if not A:
        return []
    m = len(A)
    n = len(A[0])
    H = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    col = 0
    pivots = []
    for i in range(m):
        if col >= n:
            break
        while True:
            nz = [j for j in range(col, n) if H[i][j]]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda j: abs(H[i][j]))
            j0 = nz[0]
            p = H[i][j0]
            for j in nz[1:]:
                q = H[i][j] // p
                if q:
                    for k in range(m):
                        H[k][j] -= q * H[k][j0]
                    for k in range(n):
                        U[k][j] -= q * U[k][j0]
        nz = [j for j in range(col, n) if H[i][j]]
        if not nz:
            continue
        j0 = nz[0]
        if j0 != col:
            for k in range(m):
                H[k][col], H[k][j0] = H[k][j0], H[k][col]
            for k in range(n):
                U[k][col], U[k][j0] = U[k][j0], U[k][col]
        pivots.append((i, col))
        col += 1
    y = [0] * n
    used = set()
    for (i, j) in pivots:
        s = b[i] - sum(H[i][k] * y[k] for k in used)
        if s % H[i][j]:
            return None
        y[j] = s // H[i][j]
        used.add(j)
    prows = set(i for i, _ in pivots)
    for i in range(m):
        if i in prows:
            continue
        if sum(H[i][k] * y[k] for k in used) != b[i]:
            return None
    z = [sum(U[r][k] * y[k] for k in used) for r in range(n)]
    for k in range(m):
        if sum(A[k][j] * z[j] for j in range(n)) != b[k]:
            return None
    return z


def selftest():
    import random
    from snf import solve_int
    random.seed(3)
    bad = 0
    for trial in range(400):
        m = random.randint(1, 5)
        n = random.randint(1, 5)
        A = [[random.randint(-9, 9) for _ in range(n)] for _ in range(m)]
        if random.random() < 0.5:
            z0 = [random.randint(-9, 9) for _ in range(n)]
            b = [sum(A[i][j] * z0[j] for j in range(n)) for i in range(m)]
        else:
            b = [random.randint(-30, 30) for _ in range(m)]
        r1 = solve_int(A, b)
        r2 = solve_int2(A, b)
        if (r1 is None) != (r2 is None):
            bad += 1
            print('DISAGREE', A, b, r1, r2)
    print('selftest disagreements:', bad)


if __name__ == '__main__':
    selftest()
