"""Compound knobs: pairs/triples of variables whose disturbance of the atoms OUTSIDE the
confined set cancels exactly, so the combination is a legal knob even though no single
member is.  Any such direction that leaves the current knob lattice would break a congruence.
"""
import pickle, collections, itertools, sys, time
import lib as L, model as MD, opt

v0 = opt.init()
P = L.P
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])
A = set(MD.confined_atoms(S13))
mod = MD.build(S13, v0, verbose=False)
KN = set(mod['knobs'])

MAXEXTRA = int(sys.argv[1]) if len(sys.argv) > 1 else 4
C2 = MD.load_census('24')
cand = set()
for a in A:
    cand |= C2['rev'].get(a, set())
print(f'candidate variables touching A: {len(cand)}')

info = {}
for x in sorted(cand):
    if x in KN:
        continue
    v1, t1 = MD.move(v0, {x: v0[x] + 1}, A)
    d1 = {a: t1[a] - MD.BASEP[a] for a in t1}
    ex = {a: c for a, c in d1.items() if a not in A}
    if not ex or len(ex) > MAXEXTRA:
        continue
    v5, t5 = MD.move(v0, {x: v0[x] + 5}, A)
    d5 = {a: t5[a] - MD.BASEP[a] for a in t5}
    if any(d5.get(a, 0) != 5 * c for a, c in d1.items()) or set(d5) != set(d1):
        continue                      # nonlinear
    info[x] = (frozenset(ex), ex, {a: c for a, c in d1.items() if a in A})
print(f'linear near-knobs with <= {MAXEXTRA} extra atoms: {len(info)}')

groups = collections.defaultdict(list)
for x, (k, ex, inn) in info.items():
    groups[k].append(x)
print('groups with >= 2 members:')
found = 0
for k, xs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    if len(xs) < 2:
        continue
    print(f'  extra atoms {sorted(k)}: vars {xs}')
    ea = sorted(k)
    # integer kernel of the |ea| x len(xs) matrix
    Mx = [[info[x][1].get(a, 0) for x in xs] for a in ea]
    from isolve import solve_int2
    # search 2-element cancellations exactly
    for i, j in itertools.combinations(range(len(xs)), 2):
        col_i = [Mx[r][i] for r in range(len(ea))]
        col_j = [Mx[r][j] for r in range(len(ea))]
        # find (p,q) != 0 with p*col_i + q*col_j = 0
        from math import gcd
        ok = True
        pq = None
        for r in range(len(ea)):
            a_, b_ = col_i[r], col_j[r]
            if a_ == 0 and b_ == 0:
                continue
            g = gcd(a_, b_)
            cand_pq = (b_ // g, -a_ // g)
            if pq is None:
                pq = cand_pq
            elif pq != cand_pq and pq != (-cand_pq[0], -cand_pq[1]):
                ok = False
                break
        if ok and pq:
            p_, q_ = pq
            if all(p_ * col_i[r] + q_ * col_j[r] == 0 for r in range(len(ea))):
                inner = collections.defaultdict(int)
                for a, c in info[xs[i]][2].items():
                    inner[a] += p_ * c
                for a, c in info[xs[j]][2].items():
                    inner[a] += q_ * c
                inner = {a: c for a, c in inner.items() if c}
                if inner:
                    found += 1
                    dc1 = (inner.get(22230, 0) + inner.get(22231, 0)) % P
                    dc2 = (inner.get(22229, 0) + 7376877 * inner.get(35762, 0)) % (7376877 * P)
                    print(f'    COMPOUND KNOB {p_}*x_{xs[i]} + {q_}*x_{xs[j]}  -> A-motion {inner}'
                          f'  breaksC1={dc1!=0} breaksC2={dc2!=0}')
print('compound knobs found:', found)
