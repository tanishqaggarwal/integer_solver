"""S11 step 37: full divisibility diagnostic for repairing a37887.

If A1 = 0 can be had with a37887 preserved, the twelve drop to 6 failing and the
total is 6 -> 39,027.  For each variable w of a37887 that occurs linearly, report
the coefficient c and the remainder, and try a CRT adjustment through a p-handle
to make c divide.
"""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
v = list(base); v[28730] = v[9413] * P
fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'A1 = 0; a37887 = {str(av[37887])[:40]}   score {score(v)}')

vs = sorted(set(L.avars[37887]))
print(f'\n{"var":>9} {"free":>5} {"deg":>4} {"c divides rest":>16}  gcd(c,rest)/c')
usable = []
for w in vs:
    deg = max(m.count(w) for m in L.polys[37887] if w in m)
    r = T.lin_parts(37887, w, v) if hasattr(T, 'lin_parts') else None
    if r is None:
        print(f'  x_{w:<7} {"y" if w in FREE else "n":>5} {deg:>4}   nonlinear')
        continue
    c, rest = r
    if c == 0:
        print(f'  x_{w:<7} {"y" if w in FREE else "n":>5} {deg:>4}   coefficient 0')
        continue
    ok = (rest % c == 0)
    g = math.gcd(abs(c), abs(rest))
    print(f'  x_{w:<7} {"y" if w in FREE else "n":>5} {deg:>4} {str(ok):>16}  '
          f'gcd {g} of |c| {abs(c) if abs(c) < 10**12 else str(abs(c))[:12]+"..."}')
    if ok: usable.append((w, -rest // c))
print(f'\ndirectly usable repair variables: {[w for w, _ in usable]}')
best = (score(v), None)
for w, tgt in usable:
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
        print(f'  repair x_{w} via x_{u}: score {s}  a37887 '
              f'{"ZERO" if at[37887]==0 else "nz"}  nonzero {nz}')
        if s > best[0]: best = (s, u)
        if s > 39026:
            T.save(tr, os.path.join(HERE, f'FIX_{s}.json'))
            print(f'    *** BEATS THE DELIVERABLE {s}')
print(f'\nbest {best}')
