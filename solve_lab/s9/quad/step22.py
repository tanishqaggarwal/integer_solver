"""Step 22: hand-driven construction on branch ab (x_2081=0,x_24601=0) and branch b."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002


def score(v, tag, quiet=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    if not quiet:
        print(f'[{tag:34s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))[:22]}')
    return nz, ng, ff


def zerofix(v, rounds=14, freeze=(2081, 24601), verbose=False):
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
        if verbose: print(f'    zerofix: x_{best[1]} for atom {best[3]} -> nz {best[0]}')
        ripple(v, {best[1]: best[2]})
    return v


for lbl, seeds in (('ab', {2081: 0, 24601: 0}), ('b', {24601: 0})):
    print('\n' + '=' * 78)
    print('BRANCH', lbl)
    v = H.load_assignment(BEST)
    ripple(v, seeds)
    score(v, 'raw')
    # 1) the mirror free-inputs that want ZERO (same recipe as branch a)
    for u in (24548, 14623, 14853, 31339):
        ripple(v, {u: 0})
    score(v, 'mirror free-inputs -> 0')
    # 2) the complementary load pins
    for u in (6418, 12553, 22162, 30213):
        ripple(v, {u: 0})
    score(v, 'complementary pins -> 0')
    # 3) generic zerofix
    v = zerofix(v, verbose=True)
    nz, ng, ff = score(v, 'after zerofix')
    # 4) targeted: 26731 and 33929 by matching the free partners
    ripple(v, {16742: v[19083]})
    d = 6348691 * (v[8778] - v[16144])
    if d % P: ripple(v, {8778: v[16144]})
    score(v, 'x_16742/x_8778 matched')
    v = zerofix(v)
    nz, ng, ff = score(v, 'zerofix again')
    # 5) close 688 via the x_8778 / x_33462 double knob
    k = (H688 - v[19083]) % P
    vv = list(v); ripple(vv, {8778: v[8778] + k, 33462: v[33462] + k, 16742: 0})
    ripple(vv, {16742: vv[19083]})
    dd_ = 8863713 * (vv[18956] - H688)
    if dd_ % P == 0: ripple(vv, {7497: dd_ // P})
    d2 = 6348691 * (vv[8778] - vv[16144])
    if d2 % P == 0: ripple(vv, {32253: d2 // P})
    score(vv, 'close 688 (dbl knob)')
    vv = zerofix(vv); score(vv, ' + zerofix')
    # 6) close 1618 via x_22649 / x_22152
    k2 = (H1618 - vv[30454]) % P
    ripple(vv, {22649: vv[22649] + k2, 22152: vv[22152] + k2})
    d3 = vv[24468] - H1618
    if d3 % P == 0: ripple(vv, {11436: d3 // P})
    d4 = 12604395 * (vv[22649] - vv[29524])
    if d4 % P == 0: ripple(vv, {14768: d4 // P})
    nz, ng, ff = score(vv, 'close 1618 too')
    vv = zerofix(vv)
    nz, ng, ff = score(vv, 'FINAL ' + lbl)
    for a in sorted(set(nz + ng)):
        print(f'    atom {a} ({len(atom2eq.get(a,[]))} eqs): {src[a][:150]}')
    H.save_assignment(vv, f'quad/final_{lbl}.json')
