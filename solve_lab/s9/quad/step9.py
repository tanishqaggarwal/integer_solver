"""Step 9: analyse the movers - booleanness, collateral, and whether the multiplier is free."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
v0 = H.load_assignment('quad/stateA2.json')
res = pickle.load(open('quad/sens.pkl', 'rb'))
base_nz = set(nz_checks(v0))

d19s = collections.Counter(r[0] for r in res.values())
d30s = collections.Counter(r[1] for r in res.values())
print('distinct d19083 deltas:', len(d19s), [(str(k)[:14], c) for k, c in d19s.most_common(5)])
print('distinct d30454 deltas:', len(d30s), [(str(k)[:14], c) for k, c in d30s.most_common(5)])
print('movers that are boolean-constrained:', sum(1 for t in res if t in boolv), '/', len(res))

print('\n--- minimal-collateral movers ---')
rank = sorted(res.items(), key=lambda kv: (len(kv[1][3]), kv[0]))
for t, (a19, a30, lin, nb, fx) in rank[:25]:
    print(f'x_{t:<6d} bool={t in boolv} d19={"Y" if a19 else "-"} d30={"Y" if a30 else "-"} '
          f'break({len(nb)})={nb}')

print('\n--- movers of x_30454 only (a19==0) ---')
for t, (a19, a30, lin, nb, fx) in sorted(res.items(), key=lambda kv: len(kv[1][3])):
    if a19 == 0:
        print(f'x_{t:<6d} bool={t in boolv} break({len(nb)})={nb}')

print('\n--- movers of x_19083 only (a30==0), fewest breaks ---')
c = [(t, r) for t, r in res.items() if r[1] == 0]
for t, (a19, a30, lin, nb, fx) in sorted(c, key=lambda kv: len(kv[1][3]))[:12]:
    print(f'x_{t:<6d} bool={t in boolv} break({len(nb)})={nb}')

# is the multiplier free (non-boolean linear response)?
print('\n--- multiplier test ---')
D19 = next(k for k in d19s if k)
D30 = next(k for k in d30s if k)
for t in [rank[0][0]] + [t for t in res if res[t][1] == 0][:1] + [t for t in res if res[t][0] == 0][:1]:
    for k in (5, 12345):
        v = list(v0)
        ripple(v, {t: v0[t] + k})
        ok19 = (v[19083] - v0[19083]) % P == k * res[t][0] % P
        ok30 = (v[30454] - v0[30454]) % P == k * res[t][1] % P
        print(f'  x_{t} +{k}: d19 linear={ok19}  d30 linear={ok30}  bool={t in boolv}')
print('\nneeded k19 = (H688 - x_19083)/D19 mod p =', (H688 - v0[19083]) * pow(D19, -1, P) % P)
print('needed k30 = (H1618 - x_30454)/D30 mod p =', (H1618 - v0[30454]) * pow(D30, -1, P) % P)
