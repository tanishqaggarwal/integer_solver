"""Step 13: exact reachability of  x_19083 = H688 (mod p)  and  x_30454 = H1618 (mod p)."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
v0 = H.load_assignment('quad/stateA2.json')
res = pickle.load(open('quad/sens.pkl', 'rb'))


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


c1, f1 = bcone([19083]); c2, f2 = bcone([30454])
print(f'cone(x_19083) vars={len(c1)} free={len(f1)} ; cone(x_30454) vars={len(c2)} free={len(f2)}')
print(f'union of free candidates scanned: {len(f1|f2)}   movers found: {len(res)}')

D19 = collections.defaultdict(list)
D30 = collections.defaultdict(list)
for t, (a19, a30, lin, nb, fx) in res.items():
    if a19: D19[a19].append(t)
    if a30: D30[a30].append(t)
print('\nx_19083 deltas (mod p):')
for d, ts in sorted(D19.items(), key=lambda kv: -len(kv[1])):
    print(f'  delta={d}  ({len(ts)} knobs)  bools={sum(1 for t in ts if t in boolv)}  eg {ts[:6]}')
print('x_30454 deltas (mod p):')
for d, ts in sorted(D30.items(), key=lambda kv: -len(kv[1])):
    print(f'  delta={d}  ({len(ts)} knobs)  bools={sum(1 for t in ts if t in boolv)}  eg {ts[:6]}')

t19 = (H688 - v0[19083]) % P
t30 = (H1618 - v0[30454]) % P
print('\nneeded shift of x_19083 mod p:', t19)
print('needed shift of x_30454 mod p:', t30)

HUGE31670 = 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
print('\n--- integer-combination reachability with bounded multiplicities ---')
for name, target, dd in (('x_19083', t19, D19), ('x_30454', t30, D30)):
    ds = sorted(dd, key=lambda d: -len(dd[d]))
    # each boolean knob may be toggled 0/1 -> multiplicity bound = #knobs of that delta
    print(f'\n{name}: target {str(target)[:30]}...')
    for d in ds:
        nb = len(dd[d]); free_int = [t for t in dd[d] if t not in boolv]
        hit = None
        for n in range(-nb, nb + 1):
            if (n * d - target) % P == 0:
                hit = n; break
        print(f'  delta {str(d)[:24]}.. knobs={nb} nonbool={free_int} '
              f'-> n in [-{nb},{nb}] with n*delta==target : {hit}')
        if free_int:
            print(f'     (NON-BOOLEAN knob {free_int} => arbitrary n, target reachable exactly'
                  f' if that knob is unpinned; n = {target*pow(d,-1,P)%P})')
    # two-delta combos
    if len(ds) >= 2:
        best = None
        for d1 in ds:
            for d2 in ds:
                if d1 >= d2: continue
                n1max, n2max = len(dd[d1]), len(dd[d2])
                for n1 in range(0, min(n1max, 200) + 1):
                    r = (target - n1 * d1) % P
                    for n2 in range(0, min(n2max, 200) + 1):
                        if (n2 * d2 - r) % P == 0:
                            best = (d1, n1, d2, n2)
                            break
                    if best: break
                if best: break
            if best: break
        print(f'  two-delta combo: {best}')

print('\n--- x_22649 / x_22152 pin chain ---')
print('HUGE31670 == original x_30454 :', HUGE31670 == v0[30454])
print('atom 2423 :', src[2423])
print('atom 31670:', src[31670][:200])
print('x_29524 = x_24601 * x_22152 ?')
v = list(v0)
for e in (0, 1):
    vv = list(v0); ripple(vv, {24601: e})
    print(f'   x_24601={e} -> x_29524={vv[29524]}  x_22152={vv[22152]}')
