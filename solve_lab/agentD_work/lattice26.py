"""Max number of the 12 equations satisfiable at the 39,026 placement, exactly over Z."""
import json, itertools, sys



def snf(A):
    """Smith normal form: return (D, U, V) with U*A*V = D, U,V unimodular."""
    import copy
    A = [row[:] for row in A]
    m, n = len(A), len(A[0])
    U = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]
        U[i], U[j] = U[j], U[i]

    def swap_cols(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]
        for r in V:
            r[i], r[j] = r[j], r[i]

    def addrow(i, j, c):      # row i += c * row j
        for k in range(n):
            A[i][k] += c * A[j][k]
        for k in range(m):
            U[i][k] += c * U[j][k]

    def addcol(i, j, c):      # col i += c * col j
        for r in A:
            r[i] += c * r[j]
        for r in V:
            r[i] += c * r[j]

    t = 0
    while t < min(m, n):
        # find pivot
        piv = None
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] and (best is None or abs(A[i][j]) < best):
                    best = abs(A[i][j])
                    piv = (i, j)
        if piv is None:
            break
        swap_rows(t, piv[0])
        swap_cols(t, piv[1])
        while True:
            changed = False
            for i in range(t + 1, m):
                if A[i][t]:
                    addrow(i, t, -(A[i][t] // A[t][t]))
                    if A[i][t]:
                        swap_rows(t, i)
                        changed = True
            for j in range(t + 1, n):
                if A[t][j]:
                    addcol(j, t, -(A[t][j] // A[t][t]))
                    if A[t][j]:
                        swap_cols(t, j)
                        changed = True
            if not changed and all(A[i][t] == 0 for i in range(t + 1, m)) and \
               all(A[t][j] == 0 for j in range(t + 1, n)):
                break
        t += 1
    return A, U, V


def solve_int(A, b):
    """Solve A x = b over Z.  A: m x n list of lists.  Return x or None."""
    m = len(A)
    n = len(A[0])
    D, U, V = snf([r[:] for r in A])
    Ub = [sum(U[i][k] * b[k] for k in range(m)) for i in range(m)]
    y = [0] * n
    for i in range(min(m, n)):
        d = D[i][i]
        if d == 0:
            if Ub[i] != 0:
                return None
            continue
        if Ub[i] % d:
            return None
        y[i] = Ub[i] // d
    for i in range(min(m, n), m):
        if Ub[i] != 0:
            return None
    x = [sum(V[r][c] * y[c] for c in range(n)) for r in range(n)]
    # verify
    for i in range(m):
        if sum(A[i][j] * x[j] for j in range(n)) != b[i]:
            return None
    return x


d = json.load(open('opt26.json'))
Eq = d['E']
COLS = d['cols']
M = d['M']
base = [int(x) for x in d['base']]
gens = [(u, [int(x) for x in g]) for u, g in d['gens']]
print('E', len(Eq), 'cols', COLS, 'gens', [u for u, _ in gens])
G = [[g[i] for u, g in gens] for i in range(len(COLS))]      # 8 x ngen

best = (0, None)
results = {}
for size in range(len(Eq), 0, -1):
    found = None
    for S in itertools.combinations(range(len(Eq)), size):
        A = []
        b = []
        for si in S:
            row = M[si]
            # sum_j row[j]*(base[j] + sum_k G[j][k]*x_k) = 0
            A.append([sum(row[j] * G[j][k] for j in range(len(COLS))) for k in range(len(gens))])
            b.append(-sum(row[j] * base[j] for j in range(len(COLS))))
        x = solve_int(A, b)
        if x is not None:
            found = (S, x)
            break
    results[size] = found is not None
    print(f'size {size}: solvable = {found is not None}' + (f'  S={[Eq[i] for i in found[0]]} x={found[1]}' if found else ''))
    if found:
        best = (size, found)
        break
print()
print('MAX satisfiable of the 12 =', best[0], ' => failing =', len(Eq) - best[0],
      ' => score =', L_NEQ if False else 39033 - (len(Eq) - best[0]))
