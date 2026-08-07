"""Corrected tuner: knobs = the FREED HANDLES ONLY.

tunediag2 established the affine model is exact at the required scale (12/12 predictions
correct at a delta of 10^728) and that the deliverable fixes 18 baseline rows while
breaking none.  So the target lies in the span of the freed handles alone.  The earlier
tuner used ~40 knobs, which let the sparse solver pick degenerate solutions that satisfy
the targeted rows but wreck others; restricting to the verified-affine freed set fixes it.

Validated on the calibration row before any candidate is priced.
"""
import sys, os, json, time, collections, itertools, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine3 as E3
import price as PR, fscore, sparse

NEQ = PR.NEQ
vd = PR.load_deliverable()
P = PR.TunedPricer(vd)


def tune(handles, budget=200.0, nprobe=10, want=False):
    t0 = time.time()
    freed, demote = PR.closure(handles)
    if freed is None:
        return {'ok': False, 'why': 'closure cap'}
    eng = E3.Eng(demote)
    seed = {f: vd[f] for f in eng.FREE if vd[f] != 0}
    for u in freed:
        if P.v_unc[u]:
            seed[u] = P.v_unc[u]
        else:
            seed.pop(u, None)
    v0 = eng.forward(seed)
    bad0 = eng.badatoms(v0)
    FAILS = sorted(fscore.fails(bad0))
    base = NEQ - len(FAILS)

    aff, cols = PR._affine_cols(eng, v0, bad0, sorted(freed))   # KNOBS = FREED ONLY
    if not aff:
        return {'ok': True, 'score': base, 'base_score': base, 'nknobs': 0,
                'freed': freed, 'secs': time.time() - t0}

    rows, rhs = [], []
    for e in FAILS:
        cm = collections.defaultdict(int); const = 0
        for c, a in H.eqt[e][2]:
            if a < 0:
                const += c
            else:
                cm[a] += c
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
        rhs.append(-(const + sum(c * bad0[a] for a, c in cm.items() if a in bad0)))

    order = [i for i in range(len(rows)) if rows[i]]
    keep, sols = [], []
    for i in order:
        if time.time() - t0 > budget * 0.75:
            break
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial
            sols.append(s)

    best = (base, None, None)
    if sols:
        idx = sorted(set([len(sols) - 1] +
                         [round(k * (len(sols) - 1) / max(1, nprobe - 1)) for k in range(nprobe)]))
        for j in idx:
            ns = dict(seed)
            for f, d in sols[j].items():
                if d:
                    ns[f] = v0[f] + d
            try:
                v = eng.forward(ns)
                sc = fscore.score(eng.badatoms(v))
            except Exception:
                continue
            if sc > best[0]:
                best = (sc, j + 1, ns)
            if time.time() - t0 > budget:
                break
    out = {'ok': True, 'score': best[0], 'base_score': base, 'nrows_target': len(order),
           'nsol': len(sols), 'at_prefix': best[1], 'nknobs': len(aff), 'freed': freed,
           'secs': time.time() - t0}
    if want and best[2] is not None:
        out['seed'] = best[2]; out['eng'] = eng
    return out


CANDS = [
    ([642, 28730, 29854, 31864],   'CALIBRATION (deliverable) c=27994 4971.va'),
    ([10509, 20157, 32245, 33044], 'c=27634 23762.va'),
    ([9541, 19546, 25227, 31891],  'c=31049 14869.vb'),
    ([3260, 11588, 30400, 37248],  'c=30609 14692.va'),
    ([2493, 3022, 6019, 15174],    'c=34711 3154.vb'),
    ([1405, 3052, 4806, 16433],    'c=35056 3154.va'),
    ([9337, 17894, 23336, 33996],  'c=6593  25642.va'),
    ([19053, 21505, 22193, 23910], 'c=25642 37991.vb'),
    ([10074, 16399, 16800, 35694], 'c=37034 11182.va'),
    ([1768, 6389, 26662, 31362],   'c=23089 37919.vb'),
    ([6254, 7439, 21115, 38560],   'c=37220 33113.va'),
    ([1079, 15006, 15333, 32131],  'c=7383  16012.vb'),
]

print('=== TUNER CALIBRATION: must recover 39,026 from the deliverable\'s own site ===',
      flush=True)
r0 = tune(CANDS[0][0], budget=240, want=True)
print(f'  base(untuned) {r0["base_score"]} -> TUNED {r0["score"]}   '
      f'({r0["nknobs"]} knobs, {r0["nsol"]} greedy sols, best at prefix {r0["at_prefix"]}, '
      f'{r0["secs"]:.0f}s)', flush=True)
CAL_OK = r0['score'] >= 39026
print(f'  TUNER CALIBRATION {"PASSED" if CAL_OK else "FAILED"}', flush=True)
if not CAL_OK:
    print('  -> tuned scores below are NOT measurements of their sites; reporting stops here.',
          flush=True)

if CAL_OK:
    print('\n=== pricing candidates (tuned) ===', flush=True)
    res = []
    t_all = time.time()
    for hs, tag in CANDS:
        rr = tune(hs, budget=240, want=True)
        res.append((hs, tag, rr))
        print(f'  {tag:42s} {str(hs):32s} base {rr["base_score"]} -> TUNED {rr["score"]} '
              f'({rr["secs"]:.0f}s)', flush=True)
        if rr['score'] > 39026 and rr.get('eng') is not None:
            v = rr['eng'].forward(rr['seed'])
            json.dump({f"x_{i}": int(v[i]) for i in range(PR.NV) if v[i] != 0},
                      open(f'M_site_{rr["score"]}.json', 'w'))
            print(f'     *** ABOVE BASELINE {rr["score"]} -> M_site_{rr["score"]}.json ***',
                  flush=True)
    pickle.dump([(h, t, {k: v for k, v in r.items() if k not in ('eng', 'seed')})
                 for h, t, r in res], open('tune2.pkl', 'wb'))
    per = (time.time() - t_all) / max(1, len(res))
    print(f'\nthroughput {per:.0f}s/candidate -> ~{3600/per:,.0f}/hour single-core')
    print('\n=== SUMMARY ===')
    for hs, tag, rr in sorted(res, key=lambda x: -x[2]['score']):
        print(f'  {rr["score"]:6d}  {str(hs):32s} {tag}')
    print('baseline to beat: 39026')
