"""S11 step 38: a QUADRATIC repair move.

a37887 is degree 2 in every variable, so solve_lin is blind to it -- but a
quadratic a*w^2 + b*w + c = 0 has an integer root whenever the discriminant is a
perfect square.  Several variables of a37887 are FREE.  This move class has never
been in any search.
"""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P

def isqrt_exact(n):
    if n < 0: return None
    r = math.isqrt(n)
    return r if r * r == n else None

def quad_parts(a, w, v):
    """atom a = A*w^2 + B*w + C"""
    A = B = C = 0
    for m, c in L.polys[a].items():
        k = m.count(w)
        t = c
        for z in m:
            if z != w: t *= v[z]
        if k == 0: C += t
        elif k == 1: B += t
        elif k == 2: A += t
        else: return None
    return A, B, C

def quad_roots(a, w, v):
    r = quad_parts(a, w, v)
    if r is None: return []
    A, B, C = r
    if A == 0:
        if B == 0 or C % B: return []
        return [-C // B]
    D = B * B - 4 * A * C
    s = isqrt_exact(D)
    if s is None: return []
    out = []
    for num in (-B + s, -B - s):
        if num % (2 * A) == 0: out.append(num // (2 * A))
    return out

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
v = list(base); v[28730] = v[9413] * P
fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'A1 = 0; score {score(v)}; a37887 nonzero: {av[37887] != 0}')
vs = sorted(set(L.avars[37887]))
print(f'\nquadratic roots of a37887 in each variable:')
best = (score(v), None)
for w in vs:
    rts = quad_roots(37887, w, v)
    if not rts: continue
    print(f'  x_{w:<7} {"FREE" if w in FREE else "gate":<5} roots '
          f'{[str(r)[:22] for r in rts]}  atoms {sorted(L.var_atoms[w])[:5]}')
    for tgt in rts:
        cands = []
        if w in FREE: cands.append((w, tgt))
        else:
            d = definer.get(w)
            if d is not None:
                vv = list(v); vv[w] = tgt
                for u in sorted(set(L.avars[d])):
                    if u == w or u not in FREE: continue
                    nv = T.solve_lin(d, u, vv)
                    if nv is not None: cands.append((u, nv))
        for u, nv in cands:
            tr = list(v); tr[u] = nv
            fwd(tr, rounds=8)
            at = L.all_atom_values(tr)
            s = score(tr)
            nz = [a for a in range(L.NA) if at[a]]
            print(f'      set x_{u} -> score {s}  a37887 '
                  f'{"ZERO" if at[37887]==0 else "nz"}  nonzero {nz}')
            if s > best[0]:
                best = (s, u)
                if s > 39026:
                    T.save(tr, os.path.join(HERE, f'QUAD_{s}.json'))
                    print(f'        *** BEATS THE DELIVERABLE {s}  SAVED')
print(f'\nbest {best}')
