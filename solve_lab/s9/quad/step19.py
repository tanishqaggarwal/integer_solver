"""Step 19: try x_8599=1 on branch a (frees x_12186 <- x_5096) and re-run the construction."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
hits = pickle.load(open('hits8599.pkl', 'rb'))


def score(v, tag, quiet=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    if not quiet:
        print(f'[{tag:34s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))[:20]}')
    return nz, ng, ff


def zerofix(v, rounds=8):
    """repeatedly zero any nonzero check that a single free input can absorb, best-first"""
    for _ in range(rounds):
        nz = nz_checks(v) + nz_gates(v)
        if not nz: return v
        best = None
        for a in nz:
            R = resid_poly.get(a, polys[a])
            for u in sorted(set(x for m in R for x in m)):
                if u not in freeset or u in (2081, 24601): continue
                c = 0; nl = False
                for m, cc in R.items():
                    if len(m) == 1 and m[0] == u: c += cc
                    elif u in m: nl = True
                if nl or c == 0: continue
                old = v[u]; v[u] = 0; rest = evalpoly(R, v); v[u] = old
                if rest % c: continue
                val = -rest // c
                if val == v[u]: continue
                vv = list(v); ripple(vv, {u: val})
                n2 = len(nz_checks(vv)) + len(nz_gates(vv))
                if best is None or n2 < best[0]: best = (n2, u, val)
        if best is None or best[0] >= len(nz): return v
        ripple(v, {best[1]: best[2]})
    return v


v0 = H.load_assignment(BEST)
base = list(v0); ripple(base, {2081: 0})
print('reference branch a (before x_8599):')
score(base, 'branch a')

results = []
for b in hits[:12]:
    v = list(base)
    ripple(v, {b: 1 - base[b]})
    nz, ng, ff = score(v, f'+x_{b} (x_8599={v[8599]},x_38170={v[38170]})', quiet=True)
    v = zerofix(v)
    nz, ng, ff = score(v, f'x_{b}: after zerofix')
    print('      x_12186 =', 'x_5096+p*h' if v[12186] == v[5096] + P * v[33612] else 'other',
          ' x_5096 free:', 5096 not in definer, ' x_19083 = x_23758:', v[19083] == v[23758])
    results.append((len(ff), b, sorted(set(nz + ng))))
results.sort()
print('\nbest x_8599 activators on branch a:')
for f, b, nz in results[:6]:
    print(f'  x_{b}: FAIL={f} residual={nz}')
