"""S11 step 31: the compensating pair, done in FRAME 3 (x_4432 detached).

There x_4432 holds its value, so a7930 never sees the move; a22231 and a37887
both depend only on (x_19964 + x_28730), which the pair preserves; and
a22230 = x_28730 - p*x_9413 moves freely, buying eq 29125.
x_19964 is a gate output, so realise it through its own definer a1461.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
def pr(a, n=140):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
print(f'a1461 (defines x_19964): {pr(1461)}')
for w in sorted(set(L.avars[1461])):
    print(f'   x_{w:<7} {"FREE" if w in FREE else "gate a"+str(definer.get(w))}'
          f'  atoms {sorted(L.var_atoms[w])}')
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'\nframe 3 base {score(base)}')

d = base[9413] * P - base[28730]          # target: x_28730 = p*x_9413
print(f'delta = {str(d)[:30]}...')
for tag, mk in [('x_28730 alone', None), ('pair via a1461 free input', 1461)]:
    v = list(base)
    v[28730] = v[28730] + d
    if mk is not None:
        vv = list(v); vv[19964] = v[19964] - d
        done = None
        for u in sorted(set(L.avars[1461])):
            if u == 19964 or u not in FREE: continue
            nv = T.solve_lin(1461, u, vv)
            if nv is not None:
                v[u] = nv; done = u; break
        if done is None:
            print(f'  {tag}: no free input of a1461 can realise it'); continue
        print(f'  {tag}: realised through x_{done}')
    fwd(v, rounds=8)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    s = score(v)
    print(f'  {tag:<28}: a22230 {"0" if av[22230]==0 else "nz"}  '
          f'a22231 {"0" if av[22231]==0 else "nz"}  a37887 '
          f'{"0" if av[37887]==0 else "nz"}  a7930 {"0" if av[7930]==0 else "nz"}')
    print(f'      nonzero {nz}   failing {len(L.failing_eqs(av))}  score {s}')
    if s > 39026:
        T.save(v, os.path.join(HERE, f'PAIR3_{s}.json'))
        print(f'      *** BEATS THE DELIVERABLE -- saved PAIR3_{s}.json')
