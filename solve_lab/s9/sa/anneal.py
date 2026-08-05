"""Simulated annealing / stochastic search over subsets of the 1,156 boolean free
inputs, scored by ACTUAL failing equations (fast objective, exactly validated).

usage: python3 sa/anneal.py <base.json> <align:0|1> <seconds> <seed> <tag>
"""
import sys, os, json, random, time, pickle, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

BASE = sys.argv[1]
ALIGN = bool(int(sys.argv[2]))
BUDGET = float(sys.argv[3])
SEED = int(sys.argv[4])
TAG = sys.argv[5]
random.seed(SEED)

VB = lib.H.load_assignment(BASE)
BASE_NZ = lib.nz_full(VB)
BASE_FAILS = lib.true_fails(VB)
print(f'base {BASE}: {len(BASE_FAILS)} failing, {len(BASE_NZ)} nonzero atoms', flush=True)

# validate the fast objective on the base
ff = lib.fast_fails(VB, BASE_NZ)
assert sorted(ff) == sorted(BASE_FAILS), (ff, BASE_FAILS)
print('fast objective matches on base', flush=True)

K1 = lib.K1
BFREE = lib.bfree
NEVAL = 0


def decode(bits):
    v = list(VB)
    seeds = {b: 1 for b in bits if VB[b] != 1}
    ch = {}
    if ALIGN:
        seeds[5096] = K1
        seeds[33612] = 0
    if seeds:
        c, _ = lib.ripple(v, seeds)
        ch.update(c)
    if ALIGN:
        for s in ({14853: v[12186]},
                  {7068: v[2099] + 7376877 * v[642], 4432: v[19964] + v[28730]},
                  {24548: v[25442]}):
            c, _ = lib.ripple(v, s)
            ch.update(c)
    return v, set(ch)


def score(bits):
    global NEVAL
    NEVAL += 1
    v, ch = decode(bits)
    nz = lib.nz_incremental(BASE_NZ, ch, v)
    f = lib.fast_fails(v, nz)
    return len(f), f, v


# spot-check the fast objective against the real checker on random states
random.seed(SEED + 1000)
for _ in range(3):
    bs = random.sample(BFREE, random.randint(1, 3))
    n, f, v = score(bs)
    t = lib.true_fails(v)
    assert sorted(f) == sorted(t), (bs, f, t)
print('fast objective validated against harness.evaluate on 3 random states', flush=True)
random.seed(SEED)

# classify bits: "active" ones change the failing set / score at all
ts = time.time()
ACTIVE = []
for b in BFREE:
    v, ch = decode([b])
    nz = lib.nz_incremental(BASE_NZ, ch, v)
    if nz != BASE_NZ or len(lib.fast_fails(v, nz)) != len(BASE_FAILS):
        ACTIVE.append(b)
print(f'active bits: {len(ACTIVE)} of {len(BFREE)}  ({time.time()-ts:.0f}s)', flush=True)
POOL = ACTIVE if ACTIVE else BFREE

BEST = (len(BASE_FAILS), frozenset())
cur = frozenset()
curE = len(BASE_FAILS)
visited = {cur: curE}
t0 = time.time()
T0, T1 = 6.0, 0.25
restarts = 0
hits = 0

while time.time() - t0 < BUDGET:
    frac = (time.time() - t0) / BUDGET
    T = T0 * (T1 / T0) ** frac
    # neighbourhood: flip 1-3 bits
    k = random.choice([1, 1, 1, 2, 2, 3])
    nb = set(cur)
    for _ in range(k):
        b = random.choice(POOL if random.random() < 0.85 else BFREE)
        if b in nb:
            nb.discard(b)
        else:
            nb.add(b)
    nb = frozenset(nb)
    if nb in visited:
        e = visited[nb]
    else:
        e, f, v = score(nb)
        visited[nb] = e
        if e <= len(BASE_FAILS):
            hits += 1
            with open('sa/hits.jsonl', 'a') as fh:
                fh.write(json.dumps({'method': 'anneal', 'tag': TAG, 'base': BASE,
                                     'align': ALIGN, 'bits': sorted(nb),
                                     'nfail': e, 'failing': f}) + '\n')
        if e < BEST[0]:
            BEST = (e, nb)
            out = f'sa/anneal_{TAG}_{e}.json'
            lib.H.save_assignment(v, out)
            tf = lib.true_fails(v)
            print(f'*** NEW BEST {e} (checker-model {len(tf)}) bits={sorted(nb)} -> {out}',
                  flush=True)
    if e <= curE or random.random() < math.exp(-(e - curE) / T):
        cur, curE = nb, e
    if random.random() < 0.002:            # restart
        restarts += 1
        cur, curE = frozenset(), len(BASE_FAILS)

print(f'[{TAG}] evals={NEVAL} distinct={len(visited)} restarts={restarts} '
      f'states<=base={hits} best={BEST[0]} bits={sorted(BEST[1])} '
      f'time={time.time()-t0:.0f}s', flush=True)
import collections
c = collections.Counter(visited.values())
print(f'[{TAG}] score histogram (lowest 10): {sorted(c.items())[:10]}', flush=True)
