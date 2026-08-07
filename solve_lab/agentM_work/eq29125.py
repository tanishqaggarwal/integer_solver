"""Characterise the divisibility obstruction on equation 29125, exactly.

Questions to answer:
  Q1  What quantity must divide what?
  Q2  Which knobs move either side?
  Q3  Is it a property of eq 29125 alone, or of the 86-row window?
  Q4  Is the 162-knob "core infeasible" a property of the instance or of the widening?
  Q5  Are the 5 newly-freed definer vars in the knob set?

Single-row solvability is exact and window-independent:
    sum_f coef_f * d_f = -s0   is solvable over Z  <=>  gcd_f(coef_f)  divides  s0.
So if the gcd test fails for row 29125 with a given knob set, the obstruction belongs to
equation 29125 + that knob set, NOT to the choice of window.
"""
import sys, os, json, time, math, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E_broken
import engine2 as E2, fast2, mcore2 as M, chan2 as C

P = M.P
E = 29125
vd = M.load_vec()
base = E2.seed_of(vd)
v0 = E2.forward(base)
bad0 = E2.badatoms(v0)

issq, outer, terms = H.eqt[E]
print(f'=== equation {E} ===')
print(f'issq={issq}  outer={outer}  n_terms={len(terms)}')
cmap = collections.defaultdict(int)
const = 0
for c, a in terms:
    if a < 0:
        const += c
    else:
        cmap[a] += c
print(f'constant part = {const}')
print(f'{len(cmap)} distinct atoms in the equation')
print(f'atoms of eq {E} that are currently NONZERO:')
s0 = const
for a in sorted(cmap):
    if a in bad0:
        contrib = cmap[a] * bad0[a]
        s0 += contrib
        print(f'   atom {a}: coeff {cmap[a]}  value {str(bad0[a])[:36]}...')
print(f'\ns0 (current residual of eq {E}) = {s0}')
print(f'   |s0| has {len(str(abs(s0)))} digits')


def vp(n, p):
    """p-adic valuation."""
    if n == 0:
        return None
    k = 0
    while n % p == 0:
        n //= p; k += 1
    return k


print(f'   v_P(s0) = {vp(s0, P)}   (P = the 256-bit prime)')
print(f'   s0 mod P = {s0 % P}')

# ---------- knob sets ----------
FREESET = set(E2.FREE)


def knobset(tag, cand):
    cand = sorted(f for f in cand if f in FREESET)
    aff, cols = C.affine_cols(v0, bad0, cand)
    nonaff = [f for f in cand if f not in cols]
    return tag, aff, cols, nonaff, cand


sets = []

# (1) the 19-knob set of eqsolve.py
c1 = set(E2.PIN) | set(C.CLUSTERKN)
for a in sorted(bad0):
    c1 |= set(E_broken.cone(a)[1])
sets.append(knobset('eqsolve 19-knob', c1))

# (2) the 162-knob set of eqsolve2.py (2 rounds of cone widening)
c2 = set(c1)
fails0 = set(E2.eqfails(bad0))
for e in sorted(fails0):
    for c, a in H.eqt[e][2]:
        if a >= 0:
            c2 |= set(E_broken.cone(a)[1])
seen = set()
allc = set(c2)
for rnd in range(3):
    newc = [f for f in allc if f not in seen and f in FREESET]
    if not newc:
        break
    seen |= set(newc)
    aff, cols = C.affine_cols(v0, bad0, newc)
    ext = set()
    for f in aff:
        for a in cols[f]:
            ext |= set(E_broken.cone(a)[1])
    allc |= ext
    if len(seen) > 3000:
        break
sets.append(knobset('eqsolve2 widened', seen))

# (3) EVERY free var in the cone of EVERY atom of equation 29125 -- the widest set that
#     can possibly move this row at all
c3 = set()
for a in sorted(cmap):
    try:
        c3 |= set(E_broken.cone(a)[1])
    except Exception:
        pass
c3 |= set(E2.PIN)
sets.append(knobset('all cone(eq29125 atoms)', c3))

print('\n' + '=' * 78)
for tag, aff, cols, nonaff, cand in sets:
    print(f'\n--- knob set: {tag} ---')
    print(f'  candidates {len(cand)}, affine {len(aff)}, NON-affine (excluded) {len(nonaff)}')
    print(f'  Q5  the 5 freed definer vars {E2.PIN} present as affine knobs: '
          f'{ {f: (f in cols) for f in E2.PIN} }')
    coefs = {}
    for f in aff:
        co = 0
        for a, d in cols[f].items():
            c = cmap.get(a)
            if c:
                co += c * d
        if co:
            coefs[f] = co
    print(f'  Q2  knobs with a NONZERO coefficient on eq {E}: {len(coefs)}')
    if not coefs:
        print(f'      -> NO affine knob moves eq {E} at all.')
        continue
    g = 0
    for co in coefs.values():
        g = math.gcd(g, abs(co))
    print(f'  Q1  gcd of those coefficients g = {str(g)[:60]}{"..." if len(str(g))>60 else ""}')
    print(f'      |g| digits = {len(str(g))},  v_P(g) = {vp(g, P)}')
    ok = (s0 % g == 0)
    print(f'      REQUIREMENT: g | s0   ->  s0 % g = {str(s0 % g)[:50]}   SOLVABLE: {ok}')
    if not ok:
        print(f'      -> row {E} ALONE is unsatisfiable with this knob set. '
              f'Window-independent.')
    else:
        print(f'      -> row {E} alone IS satisfiable; any failure is a JOINT/window effect.')
    top = sorted(coefs.items(), key=lambda kv: -abs(kv[1]))[:8]
    print('      largest-|coef| knobs:', [(f, len(str(abs(c)))) for f, c in top],
          '(f, #digits of coef)')
    pin_c = {f: coefs.get(f) for f in E2.PIN}
    print(f'      coefficients of the 5 freed vars: '
          f'{ {f: (len(str(abs(c))) if c else 0) for f, c in pin_c.items()} } (#digits, 0 = no effect)')

print('\n' + '=' * 78)
print('NOTE: non-affine candidates are excluded from the gcd test, so each verdict is')
print('"unsatisfiable by AFFINE motion of this knob set", not "unsatisfiable outright".')
