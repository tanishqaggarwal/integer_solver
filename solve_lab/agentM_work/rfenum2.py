"""Step 2, fast: price ALL 2^8 block-subsets from a neutral base, root-firing or not.

RAW score for every subset (exact re-propagation via resid_delta), then the
simultaneous repair on only the most promising ones.  Raw is the meaningful
instrument here: the 39,026 deliverable is a TUNED-handle point, and simsolve
replaces the tuned handles, so it cannot represent that kind of optimum.
"""
import sys, os, json, time, pickle, itertools, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2, mcore2 as M, chan2 as C, fscore

NEQ = 39033
blocks = json.load(open('blocks8.json'))['blocks']
vd = M.load_vec()
seed_d = E2.seed_of(vd)
DELIV_ON = [f for f in M.bools() if vd[f] != 0]
base = dict(seed_d)
for f in DELIV_ON:
    base.pop(f, None)

v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
print(f'NEUTRAL BASE score {fscore.score(bad0)}', flush=True)
print(f'deliverable ON {DELIV_ON}: 24601 in block0 (A-slot,178), 2081 in block2 (B-slot,21)', flush=True)

NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
reps = [[f for f in b if v0[f] == 0][:NREP] for b in blocks]
print('reps:', reps, flush=True)

results = {}
t0 = time.time()
for mask in itertools.product([0, 1], repeat=8):
    choices = [reps[i] if mask[i] else [None] for i in range(8)]
    for combo in itertools.product(*choices):
        on = tuple(f for f in combo if f is not None)
        ch = {f: 1 for f in on}
        if ch:
            b1, _ = fast2.resid_delta(v0, bad0, ch)
        else:
            b1 = bad0
        raw = fscore.score(b1)
        results[on] = (sum(mask), bool(mask[0]) and any(mask[1:]), raw, len(b1))
print(f'{len(results)} configurations, {time.time()-t0:.0f}s', flush=True)
pickle.dump(results, open('rfenum2.pkl', 'wb'))

print('\n=== RAW score vs number of live blocks ===')
by_n = collections.defaultdict(list)
for on, (n, rf, raw, nb) in results.items():
    by_n[n].append(raw)
print(' n_live | best  | worst | mean')
for n in sorted(by_n):
    v = by_n[n]
    print(f'  {n:2d}    | {max(v)} | {min(v)} | {sum(v)//len(v)}')

print('\n=== root-firing vs 78-side-only (RAW) ===')
for flag in (False, True):
    v = [raw for on, (n, rf, raw, nb) in results.items() if rf == flag]
    print(f'  rootfire={int(flag)}: best {max(v)}, worst {min(v)}, n={len(v)}')

# restricted to exactly 2 live blocks, split by whether block0 is one of them
print('\n=== exactly 2 live blocks ===')
for flag in (False, True):
    v = [raw for on, (n, rf, raw, nb) in results.items() if n == 2 and rf == flag]
    if v:
        print(f'  rootfire={int(flag)}: best {max(v)}, worst {min(v)}, n={len(v)}')

order = sorted(results.items(), key=lambda kv: -kv[1][2])
print('\ntop 15 raw:', [(on, r[2]) for on, r in order[:15]])
print('DELIVERABLE (tuned handles, same 2 slots): 39026')

# simsolve on the top few
print('\n=== simultaneous repair on top 8 raw configs ===', flush=True)
best = (0, None)
for on, r in order[:6]:
    s = dict(base)
    for f in on:
        s[f] = 1
    t1 = time.time()
    try:
        res = C.simsolve(s)
        sc = NEQ - res[0] if res else None
    except Exception as e:
        sc = f'ERR {type(e).__name__}'
    print(f'  on={on} raw={r[2]} repaired={sc} ({time.time()-t1:.0f}s)', flush=True)
    if isinstance(sc, int) and sc > best[0]:
        best = (sc, on)
        if sc > 39026:
            C.dump(E2.forward(res[1]), f'M_rf_{sc}.json')
            print(f'  *** ABOVE BASELINE {sc} ***', flush=True)
print('\nbest repaired:', best, ' baseline 39026', flush=True)
