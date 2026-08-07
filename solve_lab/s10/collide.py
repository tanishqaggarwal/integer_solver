"""S11 step 50: do the 21 collisions admit a consistent scale?

The rational kernel is 1-dimensional: atom values are lambda*z.  Where one
variable u owns two atoms a,b, setting u by delta gives a = d_a*delta and
b = d_b*delta, so we need  z_a*d_b == z_b*d_a  EXACTLY -- lambda cancels.  Either
every collision satisfies it (and the whole 69-atom vector is realisable, giving a
FULL SOLUTION) or the first failure is the precise obstruction.
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
FREE = set(ad.FREE)
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SUPP]))
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
ps = set(piv); fc = [j for j in range(m) if j not in ps][0]
z = [Fraction(0)] * m; z[fc] = Fraction(1)
for i, pj in enumerate(piv): z[pj] = -M[i][fc]
den = 1
for x in z: den = den * x.denominator // math.gcd(den, x.denominator)
zi = [int(x * den) for x in z]
g = 0
for x in zi: g = math.gcd(g, abs(x))
zi = [x // g for x in zi]
Z = {SUPP[j]: zi[j] for j in range(m)}

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
owner = collections.defaultdict(list)
for a in SUPP:
    got = None
    for w in sorted(set(L.avars[a])):
        if w in FREE and dz(a, w): got = w; break
    if got is None:
        for w in sorted(set(L.avars[a])):
            d0 = L.definer.get(w)
            if d0 is None: continue
            for u in sorted(set(L.avars[d0])):
                if u in FREE and dz(d0, u) and dz(d0, u) % P == 0: got = u; break
            if got: break
    if got is not None: owner[got].append(a)
coll = {u: l for u, l in owner.items() if len(l) > 1}
print(f'{len(coll)} colliding variables; checking z_a*d_b == z_b*d_a for each')
ok = bad = 0
fails = []
for u, lst in sorted(coll.items()):
    a, b = lst[0], lst[1]
    # exact derivative of each atom w.r.t. u, following through the definer chain
    def eff(atom, var):
        d = dz(atom, var)
        if d: return d
        for w in set(L.avars[atom]):
            d0 = L.definer.get(w)
            if d0 is None: continue
            if var in set(L.avars[d0]):
                dw = dz(d0, var); dt = dz(d0, w)
                if dt: return dz(atom, w) * (-dw) // dt if (dz(atom, w) * dw) % dt == 0 else None
        return None
    da, db = eff(a, u), eff(b, u)
    if da is None or db is None:
        fails.append((u, a, b, 'derivative not exact')); bad += 1; continue
    lhs, rhs = Z[a] * db, Z[b] * da
    if lhs == rhs: ok += 1
    else:
        bad += 1
        fails.append((u, a, b, f'z_{a}*d_{b} != z_{b}*d_{a}'))
print(f'  consistent: {ok}   INCONSISTENT: {bad}')
for f in fails[:8]: print(f'    x_{f[0]}: atoms {f[1]},{f[2]} -> {f[3]}')
if bad == 0:
    print('\n*** ALL COLLISIONS CONSISTENT -- the 69-atom vector is realisable')
