"""2^16 re-priced with the COFACTOR knobs added — an axis the enumeration never touched.

What the enumeration has been doing: a subset W names the p-handles whose relation h = p*u is
broken; `ieng.site` frees exactly those handles plus their definer-level collateral, and
`ieng.tune` then solves in equation space over THOSE variables only.  For the witness that is
five knobs: [642, 7068, 28730, 29854, 31864].

What it has NOT been doing: the 12 cofactors are ALL free inputs (checked: every one of
[105,1329,3387,5081,5676,9413,10903,11436,14393,14768,17325,22820] is in H.FREE), so they are
never in the closure and never became knobs.  T's calibration moves them (zeroing all 12 gives
39,021 / 12 failing), and the campaign records that the cofactor freedom is 4-dimensional --
only x1329, x9413, x10903, x17325 move anything.  So this is an orthogonal direction that every
number I have reported so far holds fixed.

This script re-prices the SAME 2^16 lattice with the knob set widened to
closure(W) u COFACTORS.  Affinity is not assumed: `ieng.affine_cols` tests each added knob by
second differences and drops it if it fails, exactly as for the handles.

Usage: enumcof.py <which> <budget_secs> <nprobe> <tbud> [cofmode]
  cofmode 4  -> only the 4 cofactors that move anything (default)
  cofmode 12 -> all 12
"""
import sys, os, json, time, itertools, collections, pickle

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, sparse                                    # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else '16'
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
NPROBE = int(sys.argv[3]) if len(sys.argv) > 3 else 80
TBUD = float(sys.argv[4]) if len(sys.argv) > 4 else 180.0
COFMODE = sys.argv[5] if len(sys.argv) > 5 else '4'

COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
COF4 = [1329, 9413, 10903, 17325]
EXTRA = COF4 if COFMODE == '4' else COF12

V_UNC, BAD_UNC, CM, FAILS_UNC = ieng.V_UNC, ieng.BAD_UNC, ieng.CM, ieng.FAILS_UNC
BASE = ieng.NEQ - len(FAILS_UNC)


def tune_ext(handles, nprobe=NPROBE, budget=TBUD, want=False):
    t0 = time.time()
    freed, pin = ieng.site(handles)
    if freed is None:
        return None
    knobs = sorted(set(freed) | set(EXTRA))
    aff, cols = ieng.affine_cols(pin, knobs)
    if not aff:
        return {'score': BASE, 'changes': None, 'pin': pin, 'naff': 0}
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
        return {'score': BASE, 'changes': None, 'pin': pin, 'naff': len(aff)}
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
    return {'score': best[0], 'changes': best[1], 'pin': pin, 'naff': len(aff)}


PF = json.load(open('pfamily.json'))
SETS = {k: sorted({v['h'] for v in PF[f'incident_{k2}'].values()})
        for k, k2 in (('12', '7'), ('16', '12'), ('18', '25'))}
HL = SETS[WHICH]
NTOT = 2 ** len(HL)
D4 = [642, 28730, 29854, 31864]
CKPT = f'enumcof{WHICH}_c{COFMODE}_p{NPROBE}.pkl'

print(f'H{WHICH} = {HL}', flush=True)
print(f'  extra cofactor knobs ({COFMODE}): {EXTRA}', flush=True)
print(f'  total 2^{len(HL)} = {NTOT:,}  granularity nprobe={NPROBE} budget={TBUD}', flush=True)

cal = tune_ext(D4)
print(f'CALIBRATION on the witness: {BASE} -> {cal["score"]}  ({cal["naff"]} affine knobs)  '
      f'{"PASSED" if cal["score"] >= 39026 else "FAILED"}', flush=True)
if cal['score'] < 39026:
    print('ABORT: widening the knob set must not lose the witness.', flush=True)
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
                 'i': st['i'], 'order_len': NTOT, 'handles': HL, 'extra': EXTRA,
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
        r = tune_ext(list(W)) if W else {'score': BASE, 'changes': None}
        if r is None:
            st['errors'].append((W, 'closure cap'))
        else:
            sc = r['score']
            st['bysize'][len(W)][sc] += 1
            if sc > st['best'][0]:
                st['best'] = (sc, W)
            if sc > 39026:
                st['above'].append((W, sc))
                if r['changes']:
                    bad, v = ieng.resid(V_UNC, BAD_UNC, r['changes'], r['pin'])
                    fn = f'M_cof{WHICH}_{sc}_{"_".join(map(str, W))}.json'[:120]
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
print(f'\n=== H{WHICH} + cofactors: {st["i"]:,}/{NTOT:,} priced  ({el:.0f}s, '
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
