"""S11 step 29: a37887 = x_4432 * Q with Q linear.  Find a cheap variable in Q."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
poly = L.polys[37887]
allq = all(4432 in m for m in poly)
print(f'every monomial of a37887 contains x_4432 : {allq}')
Q = {}
for m, c in poly.items():
    rest = list(m); rest.remove(4432)
    if len(rest) != 1:
        print(f'  unexpected monomial {m}')
        continue
    Q[rest[0]] = Q.get(rest[0], 0) + c
print(f'Q is linear in {len(Q)} variables')
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
cheap = sorted(((abs(c), w, c) for w, c in Q.items()),
               key=lambda t: (t[0], len(L.var_atoms[t[1]])))
print(f'\n{"|coef|":>7} {"var":>9} {"free":>6} {"atoms":>6}  value')
for ac, w, c in cheap[:16]:
    print(f'{c:>7} x_{w:<8} {"yes" if w in FREE else "no":>6} '
          f'{len(L.var_atoms[w]):>6}  {str(base[w])[:26]}')

# now: set A1 = 0 and repair Q through the cheapest usable variable
v = list(base); v[28730] = v[9413] * P
fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'\nafter A1=0: a37887 = {str(av[37887])[:34]}  score {score(v)}')
best = (score(v), None)
for ac, w, c in cheap:
    if w == 4432: continue
    tgt = T.solve_lin(37887, w, v)
    cands = []
    if tgt is not None and tgt != v[w]:
        if w in FREE: cands.append((w, tgt, 'direct'))
        else:
            d = definer.get(w)
            if d is not None:
                vv = list(v); vv[w] = tgt
                for u in sorted(set(L.avars[d])):
                    if u == w or u not in FREE: continue
                    nv = T.solve_lin(d, u, vv)
                    if nv is not None: cands.append((u, nv, f'via x_{w}'))
    for u, nv, how in cands:
        tr = list(v); tr[u] = nv
        fwd(tr, rounds=8)
        at = L.all_atom_values(tr)
        s = score(tr)
        nzt = [a for a in range(L.NA) if at[a]]
        print(f'  repair via x_{u:<7} ({how:<12}) coef {c:>10} -> score {s}  '
              f'a37887 {"ZERO" if at[37887]==0 else "nz"}  nonzero {nzt}')
        if s > best[0]:
            best = (s, u)
            T.save(tr, os.path.join(HERE, f'q_{s}.json'))
            if s > 39026: print(f'    *** BEATS THE DELIVERABLE: {s}')
print(f'\nbest {best}')
