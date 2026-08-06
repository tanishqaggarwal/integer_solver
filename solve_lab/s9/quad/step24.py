"""Step 24: on the x_8599 route, close 688 via the free x_21589 and attack the 2nd core."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626


def score(v, tag, quiet=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    if not quiet:
        print(f'[{tag:34s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))[:20]}')
    return nz, ng, ff


def zerofix(v, rounds=16, freeze=(2081, 24601), verbose=False):
    for _ in range(rounds):
        nz = nz_checks(v) + nz_gates(v)
        if not nz: return v
        best = None
        for a in nz:
            R = resid_poly.get(a, polys[a])
            for u in sorted(set(x for m in R for x in m)):
                if u not in freeset or u in freeze: continue
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
                if best is None or n2 < best[0]: best = (n2, u, val, a)
        if best is None or best[0] >= len(nz): return v
        if verbose: print(f'   zf: x_{best[1]} (atom {best[3]}) -> {best[0]}')
        ripple(v, {best[1]: best[2]})
    return v


v = H.load_assignment('quad/state8599.json')
score(v, 'state8599 (best x_8599 route)')
k = (H688 - v[19083]) % P
ripple(v, {21589: v[21589] + k})
print('x_19083 == H688 mod p:', (v[19083] - H688) % P == 0)
score(v, ' + x_21589 shift')
ripple(v, {16742: v[19083]})
d = 8863713 * (v[18956] - H688)
if d % P == 0: ripple(v, {7497: d // P})
score(v, ' + x_16742/handle')
v = zerofix(v, verbose=True)
nz, ng, ff = score(v, 'after zerofix')

# try the ~90 complement bits on the second core (S9 §9)
if nz:
    print('\nsecond-core repair scan (single bits):')
    base = len(ff)
    cands = [b for b in bfree if b not in (2081, 24601)]
    outs = []
    for b in cands:
        vv = list(v); ripple(vv, {b: 1 - v[b]})
        n2 = nz_checks(vv) + nz_gates(vv)
        if len(n2) < len(nz):
            live = [a for a in range(len(polys)) if evalpoly(polys[a], vv) != 0]
            f2 = H.evaluate(CODES, vv, eqs_of(live))
            outs.append((len(f2), len(n2), b, sorted(n2)))
    outs.sort()
    for f2, n2, b, lst in outs[:8]:
        print(f'   bit x_{b}: nz={n2} FAIL={f2} {lst[:12]}')
    if outs and outs[0][0] < base:
        f2, n2, b, lst = outs[0]
        ripple(v, {b: 1 - v[b]}); v = zerofix(v)
        score(v, f'applied x_{b}')
H.save_assignment(v, 'quad/state8599b.json')
