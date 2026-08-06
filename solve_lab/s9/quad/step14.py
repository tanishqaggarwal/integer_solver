"""Step 14: generic mod-p knob scanner; follow the pin cascade for x_16144 and x_22152."""
import sys, collections, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

v0 = H.load_assignment('quad/stateA2.json')
base_nz = set(nz_checks(v0))


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


def scan(targets, v0=v0, verbose=True):
    """For every free input in the union of backward cones, its mod-p delta on each target."""
    cands = set()
    for t in targets:
        cands |= bcone([t])[1]
    out = {}
    for u in sorted(cands):
        v = list(v0)
        ripple(v, {u: v0[u] + 1})
        dv = tuple((v[t] - v0[t]) % P for t in targets)
        if all(d == 0 for d in dv):
            continue
        nz = set(nz_checks(v))
        out[u] = (dv, sorted(nz - base_nz))
    if verbose:
        print(f'  targets {targets}: {len(cands)} cone-free-inputs, {len(out)} movers')
        agg = collections.defaultdict(list)
        for u, (dv, nb) in out.items():
            agg[dv].append(u)
        for dv, us in sorted(agg.items(), key=lambda kv: -len(kv[1])):
            nonb = [u for u in us if u not in boolv]
            print(f'    delta {tuple(str(d)[:20]+".." if d else "0" for d in dv)}  n={len(us)}  nonbool={nonb[:6]}')
            for u in (nonb[:3] or us[:2]):
                print(f'        x_{u} bool={u in boolv} breaks={out[u][1][:10]}')
    return out


print('### knobs for x_16144 (pin target of x_8778, atom 33929) ###')
s1 = scan([16144])
print('\n### knobs for x_22152 (pin target of x_22649, atom 2423) ###')
s2 = scan([22152])
print('\n### joint scan: x_19083, x_16144 ###')
s3 = scan([19083, 16144])
print('\n### joint scan: x_30454, x_22152 ###')
s4 = scan([30454, 22152])
pickle.dump({'16144': s1, '22152': s2, 'j19': s3, 'j30': s4}, open('quad/sens2.pkl', 'wb'))
