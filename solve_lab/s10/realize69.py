"""S11 step 49: can the 69 atom values be realised simultaneously?

The support touches only 68 equations, so over Q the 68 x 69 system has a kernel.
Get the exact integer kernel vector, find each atom's setting variable, and check
whether the handles COLLIDE.  If each atom has its own private handle, they are
independent and can all be set -- which would satisfy every equation.
"""
import os, sys, json, collections, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
J = json.load(open(os.path.join(HERE, 'kervec.json')))
SUPP = J['support']
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
FREE = set(ad.FREE)
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SUPP]))
print(f'support {len(SUPP)} atoms, {len(E)} equations')
# exact rational kernel of the 68 x 69 system
M = []
for e in E:
    mm, sq, co = L.eq_atoms[e]
    M.append([Fraction(co.get(a, 0)) for a in SUPP])
n, m = len(M), len(SUPP)
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if M[i][j] != 0), None)
    if k is None: continue
    M[r_], M[k] = M[k], M[r_]
    pv = M[r_][j]; M[r_] = [x / pv for x in M[r_]]
    for i in range(n):
        if i != r_ and M[i][j] != 0:
            f = M[i][j]
            M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
    piv.append(j); r_ += 1
print(f'rational rank {r_}, kernel dim {m - r_}')
ps = set(piv)
fcs = [j for j in range(m) if j not in ps]
best = None
for fc in fcs:
    z = [Fraction(0)] * m; z[fc] = Fraction(1)
    for i, pj in enumerate(piv): z[pj] = -M[i][fc]
    den = 1
    for x in z: den = den * x.denominator // math.gcd(den, x.denominator)
    zi = [int(x * den) for x in z]
    g = 0
    for x in zi: g = math.gcd(g, abs(x))
    if g: zi = [x // g for x in zi]
    ns = sum(1 for x in zi if x)
    seed = [SUPP[j] for j in range(m) if zi[j] and SUPP[j] in
            (22229, 22230, 35758, 35759, 35760, 35761, 35762)]
    if seed and (best is None or ns < best[0]): best = (ns, zi, seed)
print(f'integer kernel vector: support {best[0]}, seed atoms {best[2]}')
zi = best[1]
tgt = {SUPP[j]: zi[j] for j in range(m) if zi[j]}
print(f'max |coefficient| = {max(abs(x) for x in zi)}')

def dz(a, w):
    s = 0
    for mo, c in L.polys[a].items():
        k = mo.count(w)
        if k == 0: continue
        if k == 1:
            t = c
            for x in mo:
                if x != w: t *= v[x]
            s += t
        else: s += 2 * c * v[w]
    return s
print('\nsetting variable for each atom in the support:')
owner = collections.defaultdict(list)
noh = []
for a in sorted(tgt):
    got = None
    for w in sorted(set(L.avars[a])):
        if w in FREE:
            d = dz(a, w)
            if d: got = (w, d, 'direct'); break
    if got is None:
        for w in sorted(set(L.avars[a])):
            d0 = L.definer.get(w)
            if d0 is None: continue
            for u in sorted(set(L.avars[d0])):
                if u in FREE and dz(d0, u) and dz(d0, u) % P == 0:
                    got = (u, dz(d0, u), 'p-handle'); break
            if got: break
    if got is None: noh.append(a); continue
    owner[got[0]].append((a, got[2], got[1]))
coll = {u: lst for u, lst in owner.items() if len(lst) > 1}
print(f'  distinct setting variables: {len(owner)} for {len(tgt)} atoms')
print(f'  atoms with NO handle: {noh}')
print(f'  COLLIDING variables (one variable owns several atoms): {len(coll)}')
for u, lst in list(coll.items())[:10]:
    print(f'    x_{u} owns {[(a, k) for a, k, _ in lst]}')
