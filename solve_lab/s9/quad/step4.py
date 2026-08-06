"""Step 4: trace the cones of the remaining pinned vars / handles."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *
import cone

v = H.load_assignment('quad/stateA1.json')
cone.polys = polys


def pd(x):
    if x == 0: return 'ZERO'
    k = 0
    y = x
    while y % P == 0:
        y //= P; k += 1
    return (f'p^{k}*' if k else '') + f'{y.bit_length()}b'


def show(u, depth=0, maxdepth=5, seen=None):
    if seen is None: seen = set()
    pad = '  ' * depth
    a = definer.get(u)
    if u in seen:
        print(f'{pad}x_{u} (seen) = {pd(v[u])}'); return
    seen.add(u)
    if a is None:
        print(f'{pad}x_{u} = FREE{"[bool]" if u in boolv else ""}  val={pd(v[u])}'); return
    if depth > maxdepth:
        print(f'{pad}x_{u} <-[{a}] ... val={pd(v[u])}'); return
    print(f'{pad}x_{u} <-[{a}] {src[a][:110]}   val={pd(v[u])}')
    for w in sorted(set(u2 for m in polys[a] for u2 in m)):
        if w != u: show(w, depth + 1, maxdepth, seen)


for u in [14257, 18956, 32989, 24468, 7133, 6418, 25607, 12553, 15029, 22162, 36202, 30213, 19247, 15574, 38170, 32453, 23535]:
    print('=' * 70)
    show(u, maxdepth=4)
