#!/usr/bin/env python3
"""THE CANCELLATION EXPERIMENT.

min |{e : sum_a c_{e,a} v_a != 0}| over atom vectors v in the image of the atom map.

The residual cluster of the 39,026 witness has knob variables that occur ONLY in
residual atoms, so perturbing them is exactly realisable (no collateral).  This
builds the exact affine map knobs -> atom values -> equation cores and then
enumerates EVERY subset of the touched equations, testing integer solvability by
Smith normal form.  The answer is the true optimum for this support.
"""
import collections, itertools, sys, json, os
from model import Model, load_assign
import sympy
from sympy.matrices.normalforms import smith_normal_form

HERE = os.path.dirname(os.path.abspath(__file__))
P = 2**256 - 2**32 - 977
M = Model()
wit = load_assign(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
av = [M.atom_val(a, wit) for a in range(M.na)]

v2a = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        v2a[x].append(i)

# ---- residual support -------------------------------------------------------
RES = [23432, 23433, 36225, 36226, 36227, 36228, 36229, 37537, 37538]
COMP = [23434]                    # candidate free compensator
SUP = RES + COMP
print("support:", [(a, M.src[a][:60], len(M.atom_eqs[a])) for a in SUP])

# ---- knobs: variables occurring only inside the support ----------------------
allvars = set()
for a in SUP:
    allvars |= M.avars[a]
knobs = []
for x in sorted(allvars):
    outside = [a for a in v2a[x] if a not in SUP]
    if not outside:
        knobs.append(x)
print("\nfree knob variables (occur only inside the support):", knobs)

# effect of a unit change in each knob on each support atom (exact, linear:
# every atom is linear in each knob here -- verified numerically below)
D = {}
for k in knobs:
    w = list(wit)
    w[k] = wit[k] + 1
    D[k] = [M.atom_val(a, w) - av[a] for a in SUP]
    w[k] = wit[k] + 2
    d2 = [M.atom_val(a, w) - av[a] for a in SUP]
    assert all(d2[i] == 2 * D[k][i] for i in range(len(SUP))), f"knob {k} not linear"
print("all knobs act linearly on the support (checked)")

# ---- equations touched ------------------------------------------------------
E = sorted({e for a in SUP for e, _ in M.atom_eqs[a]})
print(f"\nequations touched by the support: {len(E)}  -> {E}")
outside_comp = sorted({e for e, _ in M.atom_eqs[23434]} - set(
    {e for a in RES for e, _ in M.atom_eqs[a]}))
print("equations a23434 adds beyond the 9-atom residual:", outside_comp)

# core_e(d) = base_e + sum_k Mat[e][k]*d_k
base = []
Mat = []
for e in E:
    coef = {a: c for c, a in M.eq_terms[e]}
    base.append(sum(coef.get(a, 0) * av[a] for a in SUP))
    Mat.append([sum(coef.get(SUP[i], 0) * D[k][i] for i in range(len(SUP)))
                for k in knobs])
print("\nbase cores (0 = already satisfied):")
for i, e in enumerate(E):
    print(f"   eq{e:6d}  base {'ZERO' if base[i] == 0 else 'nonzero'}   row {Mat[i]}"[:150])


def solvable(rows, rhs):
    """Does A x = b have an INTEGER solution?  Smith normal form test."""
    if not rows:
        return True
    A = sympy.Matrix(rows)
    b = sympy.Matrix(rhs)
    # augmented rank test over Q first (cheap necessary condition)
    if A.rank() != A.row_join(b).rank():
        return False
    m, n = A.shape
    # A = U S V  with U,V unimodular; solve S y = U b
    S, U, V = smith_normal_form(A, domain=sympy.ZZ), None, None
    # sympy's smith_normal_form does not return transforms in all versions;
    # fall back to a direct HNF-based solve via sympy.diophantine-free method:
    from sympy.matrices.normalforms import hermite_normal_form  # noqa
    # use the invariant-factor criterion: solvable over Z iff for every k,
    # d_k(A) divides d_k([A|b]) where d_k = gcd of k x k minors.  Use SNF of both.
    Sa = smith_normal_form(A, domain=sympy.ZZ)
    Sab = smith_normal_form(A.row_join(b), domain=sympy.ZZ)
    da = [Sa[i, i] for i in range(min(Sa.shape)) if Sa[i, i] != 0]
    dab = [Sab[i, i] for i in range(min(Sab.shape)) if Sab[i, i] != 0]
    if len(da) != len(dab):
        return False
    pa = 1; pb = 1
    for i in range(len(da)):
        pa *= da[i]; pb *= dab[i]
        if pb % pa != 0:
            return False
    return True


best = None
results = []
n = len(E)
for r in range(n, 0, -1):
    found = False
    for T in itertools.combinations(range(n), r):
        rows = [Mat[i] for i in T]
        rhs = [-base[i] for i in T]
        if solvable(rows, rhs):
            found = True
            best = (r, [E[i] for i in T])
            break
    if found:
        break
print(f"\n*** MAX simultaneously satisfiable equations in the support: {best[0]} of {n}")
print(f"    satisfied set: {best[1]}")
print(f"    => failing = {n - best[0]}  => SCORE = {M.ne - (n - best[0])}")
json.dump({'E': E, 'best_r': best[0], 'best_set': best[1],
           'score': M.ne - (n - best[0])},
          open(os.path.join(HERE, 'cancel_result.json'), 'w'))
