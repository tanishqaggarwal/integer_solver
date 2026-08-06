"""Step 2: dissect every nonzero residual on a quadrant branch."""
import sys, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

BRANCH = {'a': {2081: 0}, 'b': {24601: 0}, 'ab': {2081: 0, 24601: 0}}[sys.argv[1] if len(sys.argv) > 1 else 'a']

v0 = H.load_assignment(BEST)
v = list(v0)
ripple(v, BRANCH)
nz = nz_checks(v)
# gate atoms that are still nonzero must be counted too
ngate = nz_gates(v)
print('branch', BRANCH, 'nz checks', len(nz), 'nz gates', ngate)


def vinfo(u):
    a = definer.get(u)
    if a is None:
        return f'FREE{"(bool)" if u in boolv else ""}'
    return f'<-[{a}]'


def pdiv(x):
    if x == 0: return 'ZERO'
    k = 0
    while x % P == 0:
        x //= P; k += 1
    return f'p^{k}' + (f' * {x.bit_length()}b' if k else f'({x.bit_length()}b)')


print('\n### per-atom dissection ###')
for a in nz + ngate:
    R = resid_poly[a] if a in resid_poly else polys[a]
    r = evalpoly(R, v)
    vs = sorted(set(u for m in R for u in m))
    print(f'\n--- atom {a}  resid={pdiv(r)}  ({len(vs)} vars)')
    print(f'    src: {src[a][:260]}')
    # per-variable: can it linearly absorb?
    for u in vs:
        nvv = solve_for_r(R, u, v) if False else None
    # direct linear absorbers
    absorb = []
    for u in vs:
        c = 0; nl = False
        for m, cc in R.items():
            if len(m) == 1 and m[0] == u: c += cc
            elif u in m: nl = True
        if nl or c == 0: continue
        old = v[u]; v[u] = 0; rest = evalpoly(R, v); v[u] = old
        if rest % c: continue
        absorb.append((u, -rest // c, c))
    for u, nvv, c in absorb[:14]:
        print(f'      x_{u:<6d} {vinfo(u):<12s} coeff={c:<10d} needs {pdiv(nvv)}  cur={pdiv(v[u])}')
    if not absorb:
        print('      (no direct linear absorber)')
