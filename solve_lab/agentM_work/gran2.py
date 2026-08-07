"""RULE 9 on a SECOND granularity axis: the greedy ROW ORDER inside ieng.tune.

`ieng.tune` builds one equation row per baseline failure, then greedily extends a kept set
in a FIXED index order, keeping row i iff the system stays solvable over Z.  Different
orders reach different maximal solvable subsets, hence different solutions, hence different
scores.  So "nothing above 39,026" at the shipped order is a statement about that order.

The first granularity axis (nprobe/budget) is SATURATED and this script records why:
`idx` in tune() indexes into `sols`, whose length is at most the number of target rows
(<= |FAILS_UNC| = 25).  With nprobe=80 the index set already covers every element of
`sols`, so nprobe > 80 cannot add a probe.  p10 -> p80 was a real refinement; p80 -> p400
is a no-op.  Row order is therefore the axis that is still open.

Usage: gran2.py <which> <nperm> <nsample>
"""
import sys, os, json, time, random, itertools, collections, pickle

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, sparse                                    # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else '16'
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 8
NSAMP = int(sys.argv[3]) if len(sys.argv) > 3 else 400
SEED = 20260807

V_UNC, BAD_UNC, CM, FAILS_UNC = ieng.V_UNC, ieng.BAD_UNC, ieng.CM, ieng.FAILS_UNC
BASE = ieng.NEQ - len(FAILS_UNC)


def tune_perm(handles, perm, nprobe=80, budget=180.0, want=False):
    """ieng.tune with the greedy row order replaced by `perm` (a permutation of `order`)."""
    t0 = time.time()
    freed, pin = ieng.site(handles)
    if freed is None:
        return None
    aff, cols = ieng.affine_cols(pin, freed)
    if not aff:
        return {'score': BASE, 'changes': None, 'pin': pin}
    rows, rhs = [], []
    for e in FAILS_UNC:
        cm, const = CM[e]
        row = {}
        for f in aff:
            co = 0
            for a, d in cols[f].items():
                c = cm.get(a)
                if c:
                    co += c * d
            if co:
                row[f] = co
        rows.append(row)
        rhs.append(-(const + sum(c * BAD_UNC[a] for a, c in cm.items() if a in BAD_UNC)))
    order = [i for i in range(len(rows)) if rows[i]]
    if not order:
        return {'score': BASE, 'changes': None, 'pin': pin}
    if perm:
        rank = {v: i for i, v in enumerate(perm)}
        order = sorted(order, key=lambda i: rank[i])
    keep, sols = [], []
    for i in order:
        if time.time() - t0 > budget * 0.7:
            break
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial; sols.append(s)
    best = (BASE, None)
    idx = sorted(set([len(sols) - 1] +
                     [round(k * (len(sols) - 1) / max(1, nprobe - 1)) for k in range(nprobe)])) \
        if sols else []
    for j in idx:
        ch = {f: V_UNC[f] + d for f, d in sols[j].items() if d}
        if not ch:
            continue
        try:
            bad, _ = ieng.resid(V_UNC, BAD_UNC, ch, pin)
            sc = fscore.score(bad)
        except Exception:
            continue
        if sc > best[0]:
            best = (sc, ch)
        if time.time() - t0 > budget:
            break
    return {'score': best[0], 'changes': best[1], 'pin': pin}


PF = json.load(open('pfamily.json'))
SETS = {k: sorted({v['h'] for v in PF[f'incident_{k2}'].values()})
        for k, k2 in (('12', '7'), ('16', '12'), ('18', '25'))}
HL = SETS[WHICH]
D4 = [642, 28730, 29854, 31864]

# calibration: the witness must still reach 39,026 under the identity order
cal = tune_perm(D4, None)
print(f'CALIBRATION (identity order): {BASE} -> {cal["score"]}  '
      f'{"PASSED" if cal["score"] >= 39026 else "FAILED"}', flush=True)
if cal['score'] < 39026:
    sys.exit(1)

rng = random.Random(SEED)
NPROBE_ROWS = len(FAILS_UNC)
perms = [None] + [rng.sample(range(NPROBE_ROWS), NPROBE_ROWS) for _ in range(NPERM)]

# sample: every witness-superset with |W| <= 8 (the region that reaches 39,026) + uniform
sup = [tuple(sorted(set(D4) | set(c)))
       for k in range(0, 5)
       for c in itertools.combinations([h for h in HL if h not in D4], k)]
sup = [W for W in sup if len(W) <= 8]
uni = set()
while len(uni) < NSAMP:
    k = rng.randrange(0, len(HL) + 1)
    uni.add(tuple(sorted(rng.sample(HL, k))))
sample = sorted(set(sup) | uni, key=lambda W: (len(W), W))
print(f'sample: {len(sup)} witness-supersets (|W|<=8) + {len(uni)} uniform = {len(sample)} '
      f'distinct;  {len(perms)} row orders each', flush=True)

t0 = time.time()
moved = 0
above = []
delta = collections.Counter()
best = (0, None, None)
for n, W in enumerate(sample, 1):
    scs = []
    for p in perms:
        try:
            r = tune_perm(list(W), p)
            scs.append(r['score'] if r else BASE)
        except Exception:
            scs.append(BASE)
    mx = max(scs); id0 = scs[0]
    delta[mx - id0] += 1
    if mx != id0:
        moved += 1
    if mx > best[0]:
        best = (mx, W, scs.index(mx))
    if mx > 39026:
        above.append((W, mx, scs.index(mx)))
        pi = perms[scs.index(mx)]
        r = tune_perm(list(W), pi, want=True)
        if r and r['changes']:
            bad, v = ieng.resid(V_UNC, BAD_UNC, r['changes'], r['pin'])
            fn = f'M_perm{WHICH}_{mx}_{"_".join(map(str, W))}.json'[:120]
            json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0},
                      open(fn, 'w'))
            print(f'  *** ABOVE 39026: {mx} at W={W} -> {fn} ***', flush=True)
    if n % 100 == 0:
        print(f'  [{n}/{len(sample)}] {time.time()-t0:.0f}s  moved {moved}  '
              f'above {len(above)}  best {best[0]}', flush=True)

print(f'\n=== row-order granularity, {len(sample)} subsets x {len(perms)} orders, '
      f'{time.time()-t0:.0f}s ===', flush=True)
print('max-over-orders minus identity-order score:')
for d in sorted(delta):
    print(f'  delta +{d}: {delta[d]:,}')
print(f'subsets whose score MOVED : {moved:,}')
print(f'subsets now above 39,026  : {len(above)}')
print(f'BEST {best[0]} at W={best[1]} (order index {best[2]})')
json.dump({'sample': len(sample), 'nperm': len(perms), 'moved': moved,
           'delta': {str(k): v for k, v in delta.items()},
           'above': [[list(w), s, i] for w, s, i in above],
           'best': [best[0], list(best[1]) if best[1] else None, best[2]]},
          open(f'gran2_{WHICH}.json', 'w'), indent=1)
