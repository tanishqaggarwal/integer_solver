"""Step 2: channel/representative enumeration from a ROOT-FIRING base, corrected engine.

Records RAW score (leaf ON, forward, no repair) and REPAIRED score (simsolve),
because the 39,026 deliverable is NOT a simsolve output -- simsolve replaces the
tuned handle values and loses.
"""
import sys, os, json, time, pickle, collections, itertools
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2, mcore2 as M, chan2 as C, fscore

NEQ = 39033
vd = M.load_vec()
base = E2.seed_of(vd)
v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
base_score = fscore.score(bad0)
print(f'ROOT-FIRING BASE (deliverable): score {base_score}, bad {sorted(bad0)}', flush=True)

B = M.bools()
ON = [f for f in B if v0[f] != 0]
print('ON leaves:', ON, flush=True)

mode = sys.argv[1] if len(sys.argv) > 1 else 'raw'

if mode == 'raw':
    # ---- raw single-leaf scan from the root-firing base ----
    res = {}
    t0 = time.time()
    for n, f in enumerate(B):
        if v0[f] != 0:
            continue
        b1, _ = fast2.resid_delta(v0, bad0, {f: 1})
        sc = fscore.score(b1)
        res[f] = (sc, len(b1))
        if n % 64 == 0:
            print(f'  [{n}/256] {time.time()-t0:.0f}s', flush=True)
    pickle.dump(res, open('enum2_raw.pkl', 'wb'))
    order = sorted(res.items(), key=lambda kv: -kv[1][0])
    print(f'\nRAW single-leaf ON from root-firing base ({len(res)} leaves):')
    print('  best 12:', [(f, s, nb) for f, (s, nb) in order[:12]])
    print('  worst 5:', [(f, s, nb) for f, (s, nb) in order[-5:]])
    better = [(f, s) for f, (s, nb) in order if s > base_score]
    print('  ABOVE BASELINE:', better if better else 'none')
    hist = collections.Counter(s for s, nb in res.values())
    print('  score histogram:', dict(sorted(hist.items(), reverse=True)))

elif mode == 'pairs':
    # ---- raw leaf-PAIR scan, restricted to the best raw singles ----
    res = pickle.load(open('enum2_raw.pkl', 'rb'))
    top = [f for f, (s, nb) in sorted(res.items(), key=lambda kv: -kv[1][0])[:40]]
    out = {}
    t0 = time.time()
    for i, f in enumerate(top):
        for g in top[i + 1:]:
            b1, _ = fast2.resid_delta(v0, bad0, {f: 1, g: 1})
            out[(f, g)] = fscore.score(b1)
        print(f'  [{i}/{len(top)}] {time.time()-t0:.0f}s', flush=True)
    pickle.dump(out, open('enum2_pairs.pkl', 'wb'))
    order = sorted(out.items(), key=lambda kv: -kv[1])
    print('\nRAW leaf-PAIR from root-firing base, best 12:', order[:12])
    better = [(k, s) for k, s in order if s > base_score]
    print('  ABOVE BASELINE:', better[:20] if better else 'none')

elif mode == 'off':
    # ---- turning the base's OWN ON-leaves off, and single/double swaps ----
    print('\n--- effect of turning the two base ON-leaves off ---')
    for sub in [(ON[0],), (ON[1],), tuple(ON)]:
        ch = {f: 0 for f in sub}
        b1, _ = fast2.resid_delta(v0, bad0, ch)
        print(f'  off{sub}: score {NEQ-len(E2.eqfails(b1))}  nbad {len(b1)}', flush=True)

print('done', flush=True)
