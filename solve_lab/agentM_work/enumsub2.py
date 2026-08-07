"""Exhaustive enumeration over subsets of the INCIDENT p-handle set — resumable rewrite.

Same object as enumsub.py: a subset W = the p-handle atoms allowed to be nonzero, i.e. the
handles whose defining relation h = p*u is broken.  The 39,026 witness is exactly
W = {642, 28730, 29854, 31864}, |W| = 4, so the enumeration contains the known answer.

Differences from enumsub.py, all forced by the restart or by what was asked:

  * imports `shim` first, so `harness` resolves to the rebuilt model in agentM_work.
  * REAL resume.  The old checkpoint stored `n` and `last` but the script had no restart
    path; the .pkl was lost anyway.  Here the enumeration order is fixed and deterministic
    (increasing |W|, itertools.combinations within a size), the checkpoint stores the index
    into that order, and a restart replays from it.
  * per-size DISTRIBUTIONS, not just best + count@best, printed the moment a size completes.
  * errors are COUNTED and their subsets recorded, never silently skipped.  enumsub.py had a
    bare `except Exception: continue`, so a systematically failing region would have looked
    like an absent region.  A skipped subset is not a priced subset.

No ranking, no truncation, no early cutoff: per Q, an atom can be nonzero inside an equation
that still sums to zero, so incidence filters REACHABILITY, not cost, and a subset's price
cannot be bounded below by its incidence.  Every subset is priced by re-propagation.
"""
import sys, os, json, time, itertools, collections, pickle

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401  (registers harness)
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore                                            # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else '16'
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9

PF = json.load(open('pfamily.json'))
SETS = {k: sorted({v['h'] for v in PF[f'incident_{k2}'].values()})
        for k, k2 in (('12', '7'), ('16', '12'), ('18', '25'))}
HL = SETS[WHICH]
NTOT = 2 ** len(HL)
D4 = [642, 28730, 29854, 31864]
CKPT = f'enumsub{WHICH}.pkl'

print(f'H{WHICH} = {HL}', flush=True)
print(f'  witness {D4} inside: {set(D4) <= set(HL)}   total 2^{len(HL)} = {NTOT:,}', flush=True)

cal = ieng.tune(D4)
print(f'CALIBRATION on the witness: {cal["base_score"]} -> {cal["score"]}  '
      f'{"PASSED" if cal["score"] >= 39026 else "FAILED"}', flush=True)
if cal['score'] < 39026:
    sys.exit(1)

# deterministic order: increasing |W|, itertools.combinations within a size
ORDER = [W for k in range(len(HL) + 1) for W in itertools.combinations(HL, k)]
assert len(ORDER) == NTOT

st = {'bysize': collections.defaultdict(collections.Counter), 'best': (0, None),
      'above': [], 'errors': [], 'i': 0}
if os.path.exists(CKPT):
    old = pickle.load(open(CKPT, 'rb'))
    if old.get('order_len') == NTOT:
        st['bysize'] = collections.defaultdict(collections.Counter)
        for k, c in old['bysize'].items():
            st['bysize'][k] = collections.Counter(c)
        st['best'] = tuple(old['best']); st['above'] = old['above']
        st['errors'] = old['errors']; st['i'] = old['i']
        print(f'RESUMED from {CKPT} at index {st["i"]:,}/{NTOT:,}, best {st["best"][0]}',
              flush=True)

BASE = ieng.NEQ - len(ieng.FAILS_UNC)


def save(final=False):
    pickle.dump({'bysize': {k: dict(v) for k, v in st['bysize'].items()},
                 'best': st['best'], 'above': st['above'], 'errors': st['errors'],
                 'i': st['i'], 'order_len': NTOT, 'handles': HL,
                 'complete': st['i'] == NTOT and final},
                open(CKPT, 'wb'))


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
        report_size(prev_size)
        prev_size = len(W)
    try:
        r = ieng.tune(list(W)) if W else {'ok': True, 'score': BASE}
        if not r.get('ok'):
            st['errors'].append((W, r.get('why', 'not ok')))
        else:
            sc = r['score']
            st['bysize'][len(W)][sc] += 1
            if sc > st['best'][0]:
                st['best'] = (sc, W)
            if sc > 39026:
                st['above'].append((W, sc))
                r2 = ieng.tune(list(W), want=True)
                if r2.get('changes'):
                    bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r2['changes'], r2['pin'])
                    fn = f'M_sub{WHICH}_{sc}_{"_".join(map(str, W))}.json'[:120]
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
print(f'\n=== H{WHICH}: {st["i"]:,}/{NTOT:,} priced  ({el:.0f}s, '
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
