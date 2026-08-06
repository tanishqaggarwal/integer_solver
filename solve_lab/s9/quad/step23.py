"""Step 23: branch a + x_8599 activation -> 1618 via free x_5096, 688 via x_23758?"""
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
        print(f'[{tag:32s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))[:20]}')
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


def bcone(rts):
    seen = set(rts); q = collections.deque(rts); free = set()
    while q:
        u = q.popleft(); a = definer.get(u)
        if a is None: free.add(u); continue
        for w in set(x for m in polys[a] for x in m):
            if w != u and w not in seen: seen.add(w); q.append(w)
    return free


best = None
for act in (2527, 4279, 542, 1413, 2143):
    v = H.load_assignment('quad/stateA2.json')
    ripple(v, {act: 1 - v[act]})
    print(f'\n--- activator x_{act}: x_8599={v[8599]} x_38170={v[38170]} x_15298={v[15298]}')
    v = zerofix(v)
    nz, ng, ff = score(v, f'x_{act} after zerofix')
    print('    x_12186 == x_5096 + p*x_33612 :', v[12186] == v[5096] + P * v[33612],
          ' x_19083 == x_23758 :', v[19083] == v[23758])
    # close 1618 via the free x_5096
    ripple(v, {5096: H1618})
    d = v[24468] - H1618
    if d % P == 0: ripple(v, {11436: d // P})
    score(v, ' + x_5096 := H1618')
    v = zerofix(v)
    nz, ng, ff = score(v, ' + zerofix')
    if best is None or len(ff) < best[0]:
        best = (len(ff), act, sorted(set(nz + ng)), list(v))

print('\n=== best x_8599 route:', best[0], 'via x_%d' % best[1], best[2])
v = best[3]
print('knobs for x_23758 (drives x_19083 now):')
fr = bcone([23758])
r0 = v[23758]
movers = []
for u in sorted(fr):
    vv = list(v); ripple(vv, {u: v[u] + 1})
    d = (vv[23758] - r0) % P
    if d: movers.append((u, d, u in boolv))
agg = collections.Counter(d for _, d, _ in movers)
for d, n in agg.most_common():
    nb = [u for u, dd, b in movers if dd == d and not b]
    print(f'   delta {str(d)[:26]}.. n={n} nonbool={nb[:5]}')
print('   needed shift:', (H688 - v[19083]) % P)
H.save_assignment(v, 'quad/state8599.json')
