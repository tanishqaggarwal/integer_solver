"""Exhaustive enumeration with the score taken as a MAX OVER GREEDY ROW ORDERS.

`gran2.py` showed the greedy row order inside `ieng.tune` is a live granularity axis: over a
1,193-subset sample, 84% of subsets scored strictly higher under some permuted row order, by up
to +12 -- a bigger effect than the nprobe axis (max +10), which is itself saturated at nprobe=80.
So the per-subset numbers from `enumsub2.py` are LOWER BOUNDS at one fixed order.

This script makes the complete statement at the wider granularity: every subset of the lattice,
priced under NORD row orders, score = max.  Deterministic (fixed seed), resumable, per-size
distributions, errors counted.  Cost is NORD x the single-order run.

Usage: enumsub3.py <which> <budget_secs> <nord> <nprobe> <tbud>
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
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
NORD = int(sys.argv[3]) if len(sys.argv) > 3 else 3
NPROBE = int(sys.argv[4]) if len(sys.argv) > 4 else 80
TBUD = float(sys.argv[5]) if len(sys.argv) > 5 else 180.0
SEED = 20260807

V_UNC, BAD_UNC, CM, FAILS_UNC = ieng.V_UNC, ieng.BAD_UNC, ieng.CM, ieng.FAILS_UNC
BASE = ieng.NEQ - len(FAILS_UNC)
_rng = random.Random(SEED)
PERMS = [None] + [_rng.sample(range(len(FAILS_UNC)), len(FAILS_UNC)) for _ in range(NORD - 1)]


def price(handles):
    """max over row orders; returns (score, changes, pin)."""
    freed, pin = ieng.site(handles)
    if freed is None:
        return None
    aff, cols = ieng.affine_cols(pin, freed)
    if not aff:
        return (BASE, None, pin)
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
    order0 = [i for i in range(len(rows)) if rows[i]]
    if not order0:
        return (BASE, None, pin)
    best = (BASE, None)
    for perm in PERMS:
        t0 = time.time()
        if perm:
            rank = {v: i for i, v in enumerate(perm)}
            order = sorted(order0, key=lambda i: rank[i])
        else:
            order = order0
        keep, sols = [], []
        for i in order:
            if time.time() - t0 > TBUD * 0.7:
                break
            trial = keep + [i]
            s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                          verbose=False, maxcore=400, maxcorebits=5_000_000)
            if s is not None:
                keep = trial; sols.append(s)
        if not sols:
            continue
        L = len(sols) - 1
        idx = sorted(set([L] + [round(k * L / max(1, NPROBE - 1)) for k in range(NPROBE)]))
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
            if time.time() - t0 > TBUD:
                break
    return (best[0], best[1], pin)


PF = json.load(open('pfamily.json'))
SETS = {k: sorted({v['h'] for v in PF[f'incident_{k2}'].values()})
        for k, k2 in (('12', '7'), ('16', '12'), ('18', '25'))}
HL = SETS[WHICH]
NTOT = 2 ** len(HL)
D4 = [642, 28730, 29854, 31864]
CKPT = f'enumsub3_{WHICH}_o{NORD}_p{NPROBE}.pkl'

print(f'H{WHICH} = {HL}', flush=True)
print(f'  total 2^{len(HL)} = {NTOT:,}   {NORD} row orders   nprobe={NPROBE} budget={TBUD}',
      flush=True)
cal = price(D4)
print(f'CALIBRATION on the witness: {BASE} -> {cal[0]}  '
      f'{"PASSED" if cal[0] >= 39026 else "FAILED"}', flush=True)
if cal[0] < 39026:
    sys.exit(1)

ORDER = [W for k in range(len(HL) + 1) for W in itertools.combinations(HL, k)]
st = {'bysize': collections.defaultdict(collections.Counter), 'best': (0, None),
      'above': [], 'errors': [], 'i': 0}
if os.path.exists(CKPT):
    old = pickle.load(open(CKPT, 'rb'))
    if old.get('order_len') == NTOT:
        for k, c in old['bysize'].items():
            st['bysize'][k] = collections.Counter(c)
        st['best'] = tuple(old['best']); st['above'] = old['above']
        st['errors'] = old['errors']; st['i'] = old['i']
        print(f'RESUMED at index {st["i"]:,}/{NTOT:,}, best {st["best"][0]}', flush=True)


def save(final=False):
    pickle.dump({'bysize': {k: dict(v) for k, v in st['bysize'].items()},
                 'best': st['best'], 'above': st['above'], 'errors': st['errors'],
                 'i': st['i'], 'order_len': NTOT, 'handles': HL, 'nord': NORD,
                 'complete': st['i'] == NTOT and final}, open(CKPT, 'wb'))


def report_size(k):
    c = st['bysize'][k]
    n = sum(c.values()); tot = len(list(itertools.combinations(HL, k)))
    mx = max(c) if c else 0
    tag = 'COMPLETE' if n == tot else f'partial {n}/{tot}'
    print(f'  |W|={k:2d} {tag:>16}  best {mx}  count@best {c[mx]:,}', flush=True)
    print(f'         dist: ' + '  '.join(f'{s}:{c[s]:,}' for s in sorted(c, reverse=True)),
          flush=True)


t0 = time.time(); n0 = st['i']
prev_size = len(ORDER[st['i']]) if st['i'] < NTOT else None
while st['i'] < NTOT:
    W = ORDER[st['i']]
    if len(W) != prev_size:
        report_size(prev_size); prev_size = len(W)
    try:
        r = price(list(W)) if W else (BASE, None, None)
        if r is None:
            st['errors'].append((W, 'closure cap'))
        else:
            sc = r[0]
            st['bysize'][len(W)][sc] += 1
            if sc > st['best'][0]:
                st['best'] = (sc, W)
            if sc > 39026:
                st['above'].append((W, sc))
                if r[1]:
                    bad, v = ieng.resid(V_UNC, BAD_UNC, r[1], r[2])
                    fn = f'M_ord{WHICH}_{sc}_{"_".join(map(str, W))}.json'[:120]
                    json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0},
                              open(fn, 'w'))
                    print(f'  *** ABOVE 39026: {sc} at W={W} -> {fn} ***', flush=True)
    except Exception as e:
        st['errors'].append((W, repr(e)[:120]))
    st['i'] += 1
    if st['i'] % 2000 == 0:
        el = time.time() - t0
        print(f'  [{st["i"]:,}/{NTOT:,}] {el:.0f}s {(st["i"]-n0)/max(el,1e-9):.0f}/s  '
              f'best {st["best"][0]}  |W|={len(W)}  errors {len(st["errors"])}', flush=True)
        save()
    if time.time() - t0 > BUDGET:
        print(f'  [budget stop at index {st["i"]:,}]', flush=True); break

if st['i'] == NTOT:
    report_size(prev_size)
save(final=True)

el = time.time() - t0
dist = collections.Counter()
for c in st['bysize'].values():
    dist += c
print(f'\n=== H{WHICH} max-over-{NORD}-orders: {st["i"]:,}/{NTOT:,} priced  ({el:.0f}s, '
      f'{(st["i"]-n0)/max(el,1e-9):.0f}/s, complete={st["i"] == NTOT}) ===')
for k in sorted(dist, reverse=True):
    print(f'  {k}: {dist[k]:,}')
print(f'\nabove 39026: {sum(v for k, v in dist.items() if k > 39026)}')
print(f'equal 39026: {dist.get(39026, 0):,}')
print(f'errors/skips: {len(st["errors"])}')
print(f'BEST {st["best"][0]} at W={st["best"][1]}')
print('\nby support size:')
for k in sorted(st['bysize']):
    report_size(k)
