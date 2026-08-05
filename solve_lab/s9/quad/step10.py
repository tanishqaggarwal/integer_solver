"""Step 10: can the movers' own pins follow?  x_22649~x_29524 ,  x_8778~x_16144."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

v = H.load_assignment('quad/stateA2.json')


def pd(x):
    if x == 0: return 'ZERO'
    y = x; k = 0
    while y % P == 0: y //= P; k += 1
    return (f'p^{k}*' if k else '') + f'{y.bit_length()}b'


def show(u, depth=0, maxdepth=6, seen=None):
    if seen is None: seen = set()
    pad = '  ' * depth
    a = definer.get(u)
    if u in seen:
        print(f'{pad}x_{u} (seen)'); return
    seen.add(u)
    if a is None:
        print(f'{pad}x_{u} = FREE{"[bool]" if u in boolv else ""} val={pd(v[u])} occ={len(var_atoms[u])}'); return
    if depth > maxdepth:
        print(f'{pad}x_{u} <-[{a}] ...'); return
    print(f'{pad}x_{u} <-[{a}] {src[a][:100]}  val={pd(v[u])}')
    for w in sorted(set(x for m in polys[a] for x in m)):
        if w != u: show(w, depth + 1, maxdepth, seen)


for a in (2423, 33929, 26731, 688, 1618):
    print('=' * 70)
    print(f'atom {a}: {src[a][:200]}   resid={pd(evalpoly(resid_poly.get(a,polys[a]),v))}')
print()
for u in (22649, 29524, 9899, 8778, 16144, 35795, 19083, 30454, 33787):
    print('=' * 70)
    show(u, maxdepth=5)
