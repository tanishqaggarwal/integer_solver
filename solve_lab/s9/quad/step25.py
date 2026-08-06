"""Step 25: sweep all 88 x_8599 activators on branch a with the full pipeline."""
import sys, collections, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
hits = pickle.load(open('hits8599.pkl', 'rb'))


def fails(v):
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    return H.evaluate(CODES, v, eqs_of(live))


def zerofix(v, rounds=16, freeze=(2081, 24601)):
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
                if best is None or n2 < best[0]: best = (n2, u, val)
        if best is None or best[0] >= len(nz): return v
        ripple(v, {best[1]: best[2]})
    return v


base = H.load_assignment('quad/stateA2.json')
print('second core / composite equation counts:')
for a in (26733, 28438, 32342, 36185):
    print(f'   atom {a}: {len(atom2eq.get(a,[]))} eqs')

out = []
t0 = time.time()
for i, act in enumerate(hits):
    v = list(base)
    ripple(v, {act: 1 - v[act]})
    v = zerofix(v)
    ripple(v, {5096: H1618})
    d = v[24468] - H1618
    if d % P == 0: ripple(v, {11436: d // P})
    v = zerofix(v)
    k = (H688 - v[19083]) % P
    ripple(v, {21589: v[21589] + k})
    ripple(v, {16742: v[19083]})
    d = 8863713 * (v[18956] - H688)
    if d % P == 0: ripple(v, {7497: d // P})
    v = zerofix(v)
    nz = nz_checks(v) + nz_gates(v)
    ff = fails(v)
    out.append((len(ff), act, sorted(set(nz))))
    if i % 20 == 0: print(f'  {i}/{len(hits)} {time.time()-t0:.0f}s', file=sys.stderr)
out.sort()
print(f'\nsweep of {len(hits)} activators ({time.time()-t0:.0f}s) - best 10:')
for f, act, nz in out[:10]:
    print(f'   x_{act}: FAIL={f} residual({len(nz)})={nz}')
print('\nworst:', out[-1][:2])
