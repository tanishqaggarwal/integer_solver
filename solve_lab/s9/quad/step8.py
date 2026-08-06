"""Step 8: sensitivity scan - which free inputs move x_19083 and x_30454 (mod p)?"""
import sys, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
v0 = H.load_assignment('quad/stateA2.json')


def bcone(rts):
    seen = set(rts); q = collections.deque(rts); free = set()
    while q:
        u = q.popleft()
        a = definer.get(u)
        if a is None:
            free.add(u); continue
        for w in set(x for m in polys[a] for x in m):
            if w != u and w not in seen:
                seen.add(w); q.append(w)
    return seen, free


TARGETS = [19083, 30454, 16742, 12186]
c1, f1 = bcone([19083])
c2, f2 = bcone([30454])
print('cone(19083):', len(c1), 'free', len(f1))
print('cone(30454):', len(c2), 'free', len(f2))
cands = sorted(f1 | f2)
print('candidates:', len(cands))

base_nz = set(nz_checks(v0))
t0 = time.time()
res = {}
for i, t in enumerate(cands):
    for step in (1, 2):
        v = list(v0)
        ripple(v, {t: v0[t] + step})
        d19 = (v[19083] - v0[19083]) % P
        d30 = (v[30454] - v0[30454]) % P
        if step == 1:
            a19, a30 = d19, d30
            vv = v
        else:
            lin = (d19 == 2 * a19 % P) and (d30 == 2 * a30 % P)
    if a19 == 0 and a30 == 0:
        continue
    # collateral: which checks break
    nz = set(nz_checks(vv))
    res[t] = (a19, a30, lin, sorted(nz - base_nz), sorted(base_nz - nz))
    print(f'x_{t:<6d} d19083={"0" if a19==0 else str(a19)[:22]+"..":<26s} d30454={"0" if a30==0 else str(a30)[:22]+"..":<26s} '
          f'linear={lin} newbreak={res[t][3][:8]} fixed={res[t][4]}')
    if i % 200 == 0:
        print(f'   ...{i}/{len(cands)} {time.time()-t0:.0f}s', file=sys.stderr)
print(f'movers: {len(res)}  ({time.time()-t0:.0f}s)')
pickle.dump(res, open('quad/sens.pkl', 'wb'))
