"""S11 step 62: is the obstruction functional an IDENTITY?

A left-null vector y of the Jacobian means Q = sum_a y_a * a_a has zero gradient
with respect to every free input -- i.e. Q is locally constant on the manifold --
and y.r != 0 means its value there is nonzero.  If Q is EXACTLY a function of the
two boolean switches and is nonzero at all four settings, then "all atoms zero" is
impossible and the instance is INFEASIBLE.

Step 1: recover y over the 11-atom certificate support.
Step 2: form Q = sum y_a * polys[a] symbolically and see which variables survive.
"""
import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
CERT = [1436, 3576, 3578, 7930, 7932, 15456, 15462, 21617, 21619, 40065, 41507]
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
print(f'certificate atoms: {CERT}')
print(f'  values here: {[(a, "0" if av[a]==0 else "nonzero") for a in CERT]}')

G = {a: ad.grad(a, vm) for a in CERT}
cols = sorted(set().union(*[set(g) for g in G.values()]))
print(f'union of gradient supports: {len(cols)} free inputs')
# find y with sum_a y_a * grad(a) = 0   (left null of the 11 x |cols| matrix)
M = [[G[a].get(u, 0) % P for u in cols] for a in CERT]
n, m = len(M), len(cols)
A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if A[i][j]), None)
    if k is None: continue
    A[r_], A[k] = A[k], A[r_]
    inv = pow(A[r_][j], -1, P)
    A[r_] = [x * inv % P for x in A[r_]]
    for i in range(n):
        if i != r_ and A[i][j]:
            f = A[i][j]
            A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
    piv.append(j); r_ += 1
print(f'rank of the gradient matrix: {r_} of {n} atoms -> left-null dim {n - r_}')
ys = [A[i][m:] for i in range(r_, n)]
if not ys:
    print('  no left-null vector on this support'); sys.exit()
for yi, y in enumerate(ys):
    val = sum(y[i] * av[CERT[i]] for i in range(n)) % P
    print(f'  y{yi}: support {[CERT[i] for i in range(n) if y[i]]}, '
          f'Q value here = {"NONZERO" if val else "zero"}')
    if not val: continue
    # form Q = sum y_a * polys[a] symbolically, mod p
    Q = collections.defaultdict(int)
    for i, a in enumerate(CERT):
        if not y[i]: continue
        for mono, c in L.polys[a].items():
            Q[tuple(sorted(mono))] = (Q[tuple(sorted(mono))] + y[i] * c) % P
    Q = {k: c for k, c in Q.items() if c % P}
    vars_ = sorted(set().union(*[set(k) for k in Q])) if Q else []
    print(f'     Q has {len(Q)} monomials over {len(vars_)} variables')
    print(f'     variables: {vars_[:24]}{" ..." if len(vars_) > 24 else ""}')
    if set(vars_) <= {2081, 4287}:
        print('     *** Q DEPENDS ONLY ON THE TWO SWITCHES -- evaluating all four:')
        for b1 in (0, 1):
            for b2 in (0, 1):
                s = 0
                for mono, c in Q.items():
                    t = c
                    for z in mono: t = t * (b1 if z == 2081 else b2) % P
                    s = (s + t) % P
                print(f'       (x_2081,x_4287)=({b1},{b2}) -> Q = '
                      f'{"NONZERO" if s else "ZERO"}')
