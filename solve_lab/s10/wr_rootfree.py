"""WR step 17: the root-free deformation -- drop identity row 37257 (the only row
whose wire content is the bare pin alone) and get a 4-dimensional space in which
the root moves.  Identity cost 1.  Then kill the three MIXED single-equation
checks (a39084, a39417, a41278) inside that space and measure what is left."""
import os, sys, collections, json, math, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
ROOTJ = widx[26064]

DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); FW.fwd(b2)


def int_kernel(mat, n):
    m = len(mat)
    A = [r[:] for r in mat]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda c: abs(A[r][c]))
            p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m):
                        A[i][c] -= q * A[i][p0]
                    for i in range(n):
                        U[i][c] -= q * U[i][p0]
        nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
        if nz:
            piv.append(nz[0])
    return [[U[i][c] for i in range(n)] for c in range(n) if c not in piv]


def lin_form(a):
    """linear form of atom a in d (wire coordinates only); non-wire vars sit at
    their pinned values so contribute nothing."""
    f = collections.defaultdict(int)
    for m, c in L.polys[a].items():
        ws = [u for u in m if u in widx]
        if len(m) == 1 and ws:
            f[widx[m[0]]] += c
    return dict(f)


def broken(dv):
    return [e for e in RE
            if sum(c * dv[j] for j, c in rows[e].items() if dv[j])]


def measure(dv, tag):
    v = list(b2)
    for j, u in enumerate(WIRE):
        v[u] = P + dv[j]
    FW.fwd(v, rounds=10)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'{tag}: score={L.NEQ-len(fail)} failing={len(fail)} nonzero_atoms={len(nz)}',
          flush=True)
    print(f'   identity rows broken (model): {len(broken(dv))}')
    print(f'   failing eqs: {sorted(fail)}')
    return v, nz, fail


if __name__ == '__main__':
    RE2 = [e for e in RE if e != 37257]
    mat = [[rows[e].get(j, 0) for j in range(N)] for e in RE2]
    K = int_kernel(mat, N)
    print(f'kernel of the 218 rows (37257 dropped): dim {len(K)}')
    for b in K:
        for e in RE2:
            assert sum(c * b[j] for j, c in rows[e].items()) == 0
    print('   verified')
    print(f'   root coordinate values: {[len(str(abs(b[ROOTJ]))) for b in K]} digits')
    moving = [i for i, b in enumerate(K) if b[ROOTJ]]
    print(f'   basis vectors moving the root: {len(moving)} of {len(K)}')

    MIX = [39084, 39417, 41278]
    forms = {a: lin_form(a) for a in MIX}
    print('\nmixed single-equation checks and their forms:')
    for a in MIX:
        print(f'  a{a}: support {len(forms[a])}, eqs {sorted(L.atom2eq[a])}')

    # measure the raw basis vectors
    for i, b in enumerate(K):
        if b[ROOTJ]:
            measure(b, f'root-free basis {i}')
            break

    # now impose the three mixed checks = 0 inside the kernel span
    import itertools
    from fractions import Fraction
    M2 = [[Fraction(sum(f.get(j, 0) * b[j] for j in f)) for b in K] for f in
          (forms[a] for a in MIX)]
    # also require the root to move
    rowroot = [Fraction(b[ROOTJ]) for b in K]
    print(f'\n3 x {len(K)} system for the mixed checks; root row nonzero: '
          f'{any(rowroot)}')
    # solve M2 c = 0, rowroot . c != 0
    n = len(K)
    A = [r[:] for r in M2]
    piv, r_ = [], 0
    for col in range(n):
        s = next((i for i in range(r_, len(A)) if A[i][col]), None)
        if s is None:
            continue
        A[r_], A[s] = A[s], A[r_]
        inv = A[r_][col]
        A[r_] = [x / inv for x in A[r_]]
        for i in range(len(A)):
            if i != r_ and A[i][col]:
                f = A[i][col]
                A[i] = [A[i][j] - f * A[r_][j] for j in range(n)]
        piv.append(col); r_ += 1
    free = [c for c in range(n) if c not in piv]
    print(f'   rank {r_}, free directions {len(free)}')
    sols = []
    for fc in free:
        c = [Fraction(0)] * n
        c[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            c[pc] = -A[i][fc]
        den = 1
        for x in c:
            den = den * x.denominator // math.gcd(den, x.denominator)
        ci = [int(x * den) for x in c]
        g = 0
        for x in ci:
            g = math.gcd(g, abs(x))
        if g:
            ci = [x // g for x in ci]
        d = [sum(ci[i] * K[i][j] for i in range(n)) for j in range(N)]
        rootmove = d[ROOTJ] != 0
        print(f'   solution from free col {fc}: root moves {rootmove}, '
              f'a39084={sum(forms[39084].get(j,0)*d[j] for j in forms[39084])}, '
              f'a39417={sum(forms[39417].get(j,0)*d[j] for j in forms[39417])}, '
              f'a41278={sum(forms[41278].get(j,0)*d[j] for j in forms[41278])}')
        sols.append(d)
    for i, d in enumerate(sols[:3]):
        measure(d, f'root-free + mixed-checks-zeroed {i}')
