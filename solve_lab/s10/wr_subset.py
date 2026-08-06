"""WR step 8: PARTIAL wire deformation, exact price in the linear identity model.

d = (1-p)*1_T for T a subset of the 220 wire members.  Row e breaks iff
sum_{j in T} c_ej != 0.  Search T minimising the number of broken rows subject
to T containing at least one of the four residual handle multipliers.
"""
import os, sys, collections, random, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P

WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
MULT = [17499, 22665, 28961, 28599]
MI = [widx[u] for u in MULT]
ROOT = widx[26064]

# incidence for fast evaluation
colrows = collections.defaultdict(list)
for e in RE:
    for j, c in rows[e].items():
        colrows[j].append((e, c))
ei = {e: i for i, e in enumerate(RE)}


def cost(T):
    s = [0] * len(RE)
    for j in T:
        for e, c in colrows[j]:
            s[ei[e]] += c
    return sum(1 for x in s if x)


def cost_vec(T):
    s = [0] * len(RE)
    for j in T:
        for e, c in colrows[j]:
            s[ei[e]] += c
    return s


if __name__ == '__main__':
    print(f'rows {len(RE)}  wire {N}')
    print(f'per-member identity cost (d = e_u):')
    per = sorted(((cost({j}), WIRE[j]) for j in range(N)))
    print(f'   cheapest 10: {per[:10]}')
    print(f'   the four multipliers: '
          f'{[(u, cost({widx[u]})) for u in MULT]}')
    print(f'   root x_26064: {cost({ROOT})}')
    print(f'   ALL (uniform): {cost(set(range(N)))}')
    print(f'   ALL minus root: {cost(set(range(N)) - {ROOT})}')

    # local search: start from full set, flip members
    best = None
    rnd = random.Random(7)
    for trial in range(60):
        if trial == 0:
            Tset = set(range(N))
        elif trial == 1:
            Tset = set(MI)
        else:
            Tset = set(j for j in range(N) if rnd.random() < 0.5)
            Tset |= set(MI)
        cbest = cost(Tset)
        improved = True
        while improved:
            improved = False
            order = list(range(N))
            rnd.shuffle(order)
            for j in order:
                if j in MI:
                    continue
                T2 = set(Tset)
                if j in T2:
                    T2.discard(j)
                else:
                    T2.add(j)
                c2 = cost(T2)
                if c2 < cbest:
                    Tset, cbest = T2, c2
                    improved = True
        if best is None or cbest < best[0]:
            best = (cbest, set(Tset))
            print(f'  trial {trial}: broken rows {cbest}  |T|={len(Tset)}  '
                  f'root in T: {ROOT in Tset}', flush=True)
    c, Tset = best
    bad = [e for e, s in zip(RE, cost_vec(Tset)) if s]
    print(f'\nBEST partial deformation: {c} identity rows broken, |T| = {len(Tset)}')
    print(f'   broken rows: {bad}')
    print(f'   root in T: {ROOT in Tset}; multipliers in T: {[u for u in MULT]}')
    json.dump({'T': sorted(WIRE[j] for j in Tset), 'cost': c, 'bad': bad},
              open(os.path.join(HERE, 'wr_subset.json'), 'w'))
