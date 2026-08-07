"""Complete integer linear-system solvability test, plus max-integrally-zeroable-rows search.

Agent H's stageB.solve_int is INCOMPLETE: after fraction-free elimination it sets every
non-pivot coordinate to 0 and demands the pivot divisions be exact.  A system can be
integrally solvable with non-pivot coordinates nonzero, or with a different particular
solution; solve_int reports None for those.  That produces FALSE NEGATIVES, i.e. it
UNDER-counts the number of integrally zeroable rows -- exactly the quantity H's final
note says is the binding one.

Here: solvability of  M_I t = -b_I  over Z  <=>  -b_I in the integer COLUMN lattice of M_I
<=> -b_I in the integer ROW lattice of M_I^T.  Take the row Hermite normal form of M_I^T
and reduce.  That test is complete.
"""
from flint import fmpz_mat


def _hnf_rows(cols):
    """cols: list of generator vectors (each length m).  Returns echelon HNF rows (list of lists)."""
    H = fmpz_mat([list(c) for c in cols]).hnf()
    out = []
    for r in H.tolist():
        rr = [int(x) for x in r]
        if any(rr):
            out.append(rr)
    return out


def lattice_member(H, c):
    """H: HNF rows (echelon, increasing pivots).  c: vector.  True iff c in Z-rowspan(H)."""
    c = list(c)
    m = len(c)
    piv = []
    for r in H:
        j = next((k for k in range(m) if r[k]), None)
        if j is not None:
            piv.append((j, r))
    piv.sort(key=lambda z: z[0])
    for j, r in piv:
        if c[j]:
            if c[j] % r[j]:
                return False
            f = c[j] // r[j]
            if f:
                for k in range(j, m):
                    c[k] -= f * r[k]
    return not any(c)


_Q = (1 << 61) - 1   # Mersenne prime


class ZSolver:
    """Cached integer-solvability oracle for row subsets of a fixed (M, b).

    A mod-q pass runs first: an integer solution reduces to a solution mod q, so
    mod-q inconsistency PROVES integer inconsistency.  Only survivors reach the
    (expensive, big-integer) HNF test, which decides the remaining cases exactly.
    """

    def __init__(self, M, b, n):
        self.M = M          # M[i] = list of n coefficients for row i
        self.b = b
        self.n = n
        self.cache = {}
        q = _Q
        self.Mq = [[x % q for x in row] for row in M]
        self.bq = [(-x) % q for x in b]

    def _modq_consistent(self, rows):
        """Gaussian elimination mod q on [M_rows | -b_rows]; False => no integer solution."""
        q = _Q
        n = self.n
        A = [self.Mq[i] + [self.bq[i]] for i in rows]
        r = 0
        for c in range(n):
            k = None
            for i in range(r, len(A)):
                if A[i][c]:
                    k = i
                    break
            if k is None:
                continue
            A[r], A[k] = A[k], A[r]
            inv = pow(A[r][c], q - 2, q)
            Ar = A[r] = [(x * inv) % q for x in A[r]]
            for i in range(len(A)):
                if i != r and A[i][c]:
                    f = A[i][c]
                    A[i] = [(A[i][j] - f * Ar[j]) % q for j in range(n + 1)]
            r += 1
            if r == len(A):
                break
        for i in range(r, len(A)):
            if A[i][n] and not any(A[i][j] for j in range(n)):
                return False
        return True

    def solvable(self, rows):
        key = tuple(sorted(rows))
        v = self.cache.get(key)
        if v is not None:
            return v
        if not key:
            self.cache[key] = True
            return True
        M, b = self.M, self.b
        if not self._modq_consistent(key):
            self.cache[key] = False
            return False
        # quick: any row with all-zero coefficients but nonzero b is unsolvable
        for i in key:
            if b[i] and not any(M[i]):
                self.cache[key] = False
                return False
        gens = [[M[i][j] for i in key] for j in range(self.n)]
        H = _hnf_rows(gens)
        r = lattice_member(H, [-b[i] for i in key])
        self.cache[key] = r
        return r


def max_zero_rows(M, b, n, nrows, node_cap=200000, need=None, lb=0):
    """Exact maximum number of rows integrally zeroable simultaneously.

    The family of integrally solvable row-subsets is downward closed, so this is a
    maximum-"clique"-style search: at each node keep only the rows that still extend the
    current set, and prune with |cur| + |candidates| <= best.  Exact when exhaustive.
    Returns (best_size, best_subset, exhaustive_flag, tests).
    """
    Z = ZSolver(M, b, n)
    free0 = [i for i in range(nrows) if b[i] == 0 and not any(M[i])]   # always zero
    fz = set(free0)
    rest = [i for i in range(nrows) if i not in fz]
    tests = [0]
    exhaustive = [True]
    base = len(free0)
    # `lb` seeds the incumbent: the search then proves only "opt <= lb" or returns the exact
    # opt when opt > lb.  Sound for the decision question and far cheaper.
    best = [max(base, lb), list(free0)]
    seeded = best[0] > base

    def solvable(rows):
        tests[0] += 1
        return Z.solvable(rows)

    cands0 = [i for i in rest if solvable([i])]
    if need is not None and base + len(cands0) < need:
        # sound upper bound: every row of an optimal set is individually solvable
        return base + len(cands0), None, 'bound', tests[0]

    def dfs(cur, cands):
        if tests[0] > node_cap:
            exhaustive[0] = False
            return
        if base + len(cur) + len(cands) <= best[0]:
            return
        for idx in range(len(cands)):
            if tests[0] > node_cap:
                exhaustive[0] = False
                return
            if base + len(cur) + (len(cands) - idx) <= best[0]:
                return
            i = cands[idx]
            nxt = cur + [i]
            if base + len(nxt) > best[0]:
                best[0] = base + len(nxt)
                best[1] = free0 + nxt
            newc = [j for j in cands[idx + 1:] if solvable(nxt + [j])]
            if newc:
                dfs(nxt, newc)

    dfs([], cands0)
    if seeded and best[1] == free0 and best[0] == lb:
        # never improved on the seed: opt is only known to be <= lb
        return lb, None, 'le', tests[0]
    return best[0], sorted(best[1]), exhaustive[0], tests[0]


def witness_t(M, b, n, rows):
    """Recover an actual integer t with (Mt+b)_rows = 0, or None.  Uses HNF with transform
    emulated by solving over the lattice generated by unit vectors."""
    rows = sorted(rows)
    if not rows:
        return [0] * n
    m = len(rows)
    # Build [ M_I^T | I_n ] and row-HNF it: the left block gives lattice generators, the
    # right block records the integer combination of the original generators (= columns of M_I).
    aug = [[M[i][j] for i in rows] + [1 if k == j else 0 for k in range(n)] for j in range(n)]
    H = fmpz_mat(aug).hnf().tolist()
    H = [[int(x) for x in r] for r in H]
    c = [-b[i] for i in rows]
    t = [0] * n
    piv = []
    for r in H:
        j = next((k for k in range(m) if r[k]), None)
        if j is not None:
            piv.append((j, r))
    piv.sort(key=lambda z: z[0])
    for j, r in piv:
        if c[j]:
            if c[j] % r[j]:
                return None
            f = c[j] // r[j]
            for k in range(j, m):
                c[k] -= f * r[k]
            for k in range(n):
                t[k] += f * r[m + k]
    if any(c):
        return None
    for idx, i in enumerate(rows):
        if b[i] + sum(t[j] * M[i][j] for j in range(n)) != 0:
            return None
    return t


if __name__ == '__main__':
    import random
    random.seed(7)
    # self-test: random systems, compare against brute force over a small box
    ok = 0
    for trial in range(300):
        n = random.randint(1, 3)
        m = random.randint(1, 3)
        A = [[random.randint(-4, 4) for _ in range(n)] for _ in range(m)]
        tt = [random.randint(-3, 3) for _ in range(n)]
        if random.random() < 0.5:
            bb = [-sum(A[i][j] * tt[j] for j in range(n)) for i in range(m)]   # solvable
            want = True
        else:
            bb = [random.randint(-9, 9) for _ in range(m)]
            want = None
        Z = ZSolver(A, bb, n)
        got = Z.solvable(range(m))
        if want is True:
            assert got, (A, bb)
        # brute force cross-check on a box
        brute = False
        rng = range(-8, 9)
        import itertools as it
        for cand in it.product(rng, repeat=n):
            if all(bb[i] + sum(cand[j] * A[i][j] for j in range(n)) == 0 for i in range(m)):
                brute = True
                break
        if brute:
            assert got, ('missed solvable', A, bb)
        ok += 1
        # witness recovery
        if got:
            t = witness_t(A, bb, n, range(m))
            assert t is not None, ('no witness', A, bb)
    print('zsolve self-test OK on %d random systems' % ok)
