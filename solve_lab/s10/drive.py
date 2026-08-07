"""S11 step 32: can x_8731 (zero collateral) drive x_19964?"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
def pr(a, n=130):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
for a in (1459, 1460, 1461, 8261):
    print(f'a{a}: {pr(a)}')
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av0 = L.all_atom_values(base)
print(f'\nbase: x_19964 = {str(base[19964])[:26]}  a37887 = {av0[37887]}')
for u in (8731, 9118, 20492, 19892):
    v = list(base); v[u] = v[u] + 1
    fwd(v, rounds=8)
    av = L.all_atom_values(v)
    print(f'  x_{u} +1 -> x_19964 changed by {v[19964]-base[19964]}, '
          f'x_19892 by {v[19892]-base[19892]}, x_20492 by {v[20492]-base[20492]}; '
          f'a37887 {"0" if av[37887]==0 else "nz"}  a22231 '
          f'{"0" if av[22231]==0 else "nz"}')

print('\n=== the full move: x_28730 -> p*x_9413 with x_8731 compensating ===')
d = base[9413] * P - base[28730]
# a1459 relates x_19892 and x_8731 ; find the exact coefficient
v = list(base); v[8731] = v[8731] + 1
fwd(v, rounds=8)
c19964 = v[19964] - base[19964]
print(f'd(x_19964)/d(x_8731) = {c19964}')
if c19964 and d % c19964 == 0:
    step = -d // c19964
    w = list(base)
    w[28730] = w[28730] + d
    w[8731] = w[8731] + step
    fwd(w, rounds=8)
    aw = L.all_atom_values(w)
    nz = [a for a in range(L.NA) if aw[a]]
    s = score(w)
    print(f'  x_8731 += {str(step)[:26]}')
    print(f'  a22230 {"0" if aw[22230]==0 else "nz"}  a22231 '
          f'{"0" if aw[22231]==0 else "nz"}  a37887 {"0" if aw[37887]==0 else "nz"}'
          f'  a7930 {"0" if aw[7930]==0 else "nz"}')
    print(f'  nonzero {nz}  failing {len(L.failing_eqs(aw))}  score {s}')
    if s > 39026:
        T.save(w, os.path.join(HERE, f'DRIVE_{s}.json'))
        print(f'  *** BEATS THE DELIVERABLE -- saved DRIVE_{s}.json')
else:
    print(f'  not divisible: d mod c = {d % c19964 if c19964 else "c=0"}')
