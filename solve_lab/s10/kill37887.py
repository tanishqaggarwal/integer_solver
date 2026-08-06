"""S10 step 70: extend the model to kill a37887 as well.

The 8-atom model reaches 6 of 12 in the region, but a37887 = R^2 lights up and
breaks equation 8680, restoring 7.  R is a LINEAR COMBINATION OF ATOMS:

    R = a22231 + 6*a22232 + 15*a22233 - 21*a22234 + ...

and a22232..a22235 are themselves movable:
    x_23754 by d  ->  a22232 += d, a22233 -= d       (x_21279 = 0)
    x_35619 by e  ->  a22234 += e, a22235 += e
So R can be driven to 0.  Cost: a22232..a22235 drag in equations
{6494, 8687, 22563, 35561}.  Optimise the whole enlarged region exactly.
"""
import os, sys, ast, json, math, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)

# ---- decompose a37887's root into atoms -----------------------------------
src = L.atom_src[37887]
node = ast.parse(src, mode='eval').body
assert isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult)
root = node.left
idx = {s: i for i, s in enumerate(L.atom_src)}


def unfold(n):
    out = []
    while isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
        r = n.right
        c, a = None, None
        if isinstance(r, ast.BinOp) and isinstance(r.op, ast.Mult):
            try:
                c = ast.literal_eval(r.left); a = r.right
            except Exception:
                c = None
        if c is None:
            break
        out.append((c, a)); n = n.left
    out.append((1, n)); out.reverse()
    return out


R_terms = []
for c, a in unfold(root):
    key = ast.unparse(a)
    if key in idx:
        R_terms.append((c, idx[key]))
    else:
        R_terms.append((c, None))
print('a37887 root decomposes into:')
for c, a in R_terms:
    print(f'   {c:>5} * a{a}')
assert all(a is not None for c, a in R_terms), 'root has a non-atom term'
Rval = sum(c * av[a] for c, a in R_terms)
print(f'root value now = {Rval}  (a37887 = {av[37887]}, R^2 matches: '
      f'{Rval*Rval == av[37887]})')

# ---- the enlarged model ----------------------------------------------------
ATOMS = [22229, 22230, 35758, 35759, 35760, 35761, 35762, 22231,
         22232, 22233, 22234, 22235]
E = sorted(L.eqs_of_atoms(ATOMS) | set(L.atom2eq.get(37887, ())))
print(f'\nregion: {len(ATOMS)} atoms, {len(E)} equations {E}')
M = []
for e in E:
    if e in L.atom2eq.get(37887, ()):
        continue                      # handled by the R = 0 constraint
    M.append([L.eq_atoms[e][2].get(a, 0) for a in ATOMS])
EQ_LIN = [e for e in E if e not in L.atom2eq.get(37887, ())]
print(f'linear equations in the region: {len(EQ_LIN)}; '
      f'plus R = 0 for eq {sorted(L.atom2eq.get(37887, ()))}')

D0 = v[7068] - v[2099]
K = v[4432] - v[19964]
n = len(ATOMS)
Rrow = [0] * n
for c, a in R_terms:
    if a in ATOMS:
        Rrow[ATOMS.index(a)] += c


def int_kernel(mat):
    m, nn = len(mat), len(mat[0])
    A = [r[:] for r in mat]
    U = [[1 if i == j else 0 for j in range(nn)] for i in range(nn)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(nn) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda c: abs(A[r][c])); p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(nn): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(nn) if c not in piv and A[r][c] != 0]
        if nz: piv.append(nz[0])
    return [[U[i][c] for i in range(nn)] for c in range(nn) if c not in piv]


def solve_int(rows_, rhs, nn):
    m = len(rows_)
    A = [r[:] for r in rows_]
    U = [[1 if i == j else 0 for j in range(nn)] for i in range(nn)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(nn) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda c: abs(A[r][c])); p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(nn): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(nn) if c not in piv and A[r][c] != 0]
        piv.append(nz[0] if nz else None)
    w = [0] * nn; b = list(rhs)
    for r in range(m):
        c = piv[r]
        if c is None:
            if b[r] != 0: return None
            continue
        if b[r] % A[r][c]: return None
        w[c] = b[r] // A[r][c]
        for rr in range(r + 1, m): b[rr] -= A[rr][c] * w[c]
    return [sum(U[i][c] * w[c] for c in range(nn)) for i in range(nn)]


# achievability constraints (exact):
#   A1 + 7376877*A7 == D0 (mod p)          [free k*p shift of x_7068]
#   A2 + A8         == K  (mod p)          [x_28730 alone; x_4432 fixed]
#   A9 + A10 == 0                          [x_23754 moves them oppositely]
#   A11 - A12 == 0                         [x_35619 moves them together]
#   R == 0                                 [kills a37887 / eq 8680]
print('\n=== exact optimisation, 12-atom model with R = 0 ===')
best = None
for size in range(len(EQ_LIN), 0, -1):
    hit = None
    for S in itertools.combinations(range(len(EQ_LIN)), size):
        rows = [M[i] for i in S] + [Rrow, [0]*8 + [1, 1, 0, 0], [0]*8 + [0, 0, 1, -1]]
        B = int_kernel(rows)
        if not B:
            continue
        c1 = [b[0] + 7376877 * b[6] for b in B]
        c2 = [b[1] + b[7] for b in B]
        sol = solve_int([c1 + [P, 0], c2 + [0, P]], [D0 % P, K % P], len(B) + 2)
        if sol is not None:
            y = sol[:len(B)]
            A = [sum(y[j] * B[j][i] for j in range(len(B))) for i in range(n)]
            if any(A):
                hit = (S, A); break
    print(f'  size {size:>2} of {len(EQ_LIN)}: {"SOLVABLE" if hit else "no"}', flush=True)
    if hit:
        best = (size, hit); break

if best:
    size, (S, A) = best
    failing = len(EQ_LIN) - size
    print(f'\n*** satisfied {size} of {len(EQ_LIN)} linear equations, plus eq 8680 '
          f'via R=0  ->  FAILING {failing}  SCORE {L.NEQ-failing}')
    print(f'    atoms {ATOMS}')
    print(f'    values {[str(x)[:22] for x in A]}')
    json.dump({'ATOMS': ATOMS, 'EQ': EQ_LIN, 'S': [EQ_LIN[i] for i in S],
               'A': [str(x) for x in A], 'K': str(K), 'D0': str(D0),
               'R_terms': [[c, a] for c, a in R_terms]},
              open(os.path.join(HERE, 'kill37887.json'), 'w'))
    print('    saved kill37887.json')
