"""Corrected-engine analogue of E's chanenum, enumerating over the 8 measured blocks.

Base = the deliverable's seed with its two ON leaves turned OFF (neutral start).
For each of the 2^8 subsets of the 8 blocks, turn ON one representative per block,
then run the corrected simultaneous solve.  Block 0 is the root's A-slot (178);
blocks 1..7 are the B-slot (78).  A subset containing block 0 AND at least one of
1..7 is ROOT-FIRING.  This prices root-firing vs 78-side-only under one instrument.
"""
import sys, os, json, time, pickle, itertools, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2, mcore2 as M, chan2 as C

NEQ = 39033
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
blocks = json.load(open('blocks8.json'))['blocks']
vd = M.load_vec()
seed_d = E2.seed_of(vd)
DELIV_ON = [f for f in M.bools() if vd[f] != 0]

# neutral base: deliverable handles, both leaves off
base = dict(seed_d)
for f in DELIV_ON:
    base.pop(f, None)

v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
print(f'NEUTRAL BASE: score {NEQ-len(E2.eqfails(bad0))}, {len(bad0)} bad atoms', flush=True)
print(f'deliverable ON leaves {DELIV_ON}  (24601 A-side/block0, 2081 B-side/block2)', flush=True)

reps = []
for i, b in enumerate(blocks):
    cand = [f for f in b if v0[f] == 0]
    reps.append(cand[:NREP])
print('reps per block:', reps, flush=True)

results = {}
t0 = time.time()
for mask in itertools.product([0, 1], repeat=8):
    choices = [reps[i] if mask[i] else [None] for i in range(8)]
    for combo in itertools.product(*choices):
        s = dict(base)
        on = [f for f in combo if f is not None]
        for f in on:
            s[f] = 1
        # raw score
        ch = {f: 1 for f in on}
        b1, _ = fast2.resid_delta(v0, bad0, ch) if ch else (bad0, None)
        raw = NEQ - len(E2.eqfails(b1))
        # repaired score
        try:
            r = C.simsolve(s)
            rep = NEQ - r[0] if r else None
        except Exception as e:
            rep = None
        nlive = sum(mask)
        rootfire = bool(mask[0]) and any(mask[1:])
        results[tuple(on)] = (nlive, rootfire, raw, rep)
        print(f'  mask{"".join(map(str,mask))} n={nlive} rootfire={int(rootfire)} '
              f'raw={raw} repaired={rep} on={on}', flush=True)
    if time.time() - t0 > 2400:
        print('  [time budget reached]', flush=True)
        break

pickle.dump(results, open('rfenum.pkl', 'wb'))

print('\n=== SUMMARY ===', flush=True)
by_n = collections.defaultdict(list)
for on, (nlive, rf, raw, rep) in results.items():
    by_n[nlive].append((raw, rep))
print('n_live | best_raw | best_repaired')
for n in sorted(by_n):
    raws = [a for a, b in by_n[n]]
    reps_ = [b for a, b in by_n[n] if b is not None]
    print(f'  {n:2d}   | {max(raws)}   | {max(reps_) if reps_ else "-"}')

print('\nroot-firing vs not (repaired):')
for flag in (False, True):
    vals = [rep for on, (n, rf, raw, rep) in results.items() if rf == flag and rep is not None]
    rawv = [raw for on, (n, rf, raw, rep) in results.items() if rf == flag]
    if vals:
        print(f'  rootfire={int(flag)}: best repaired {max(vals)}, best raw {max(rawv)}, n={len(vals)}')

allbest = max(((rep, on) for on, (n, rf, raw, rep) in results.items() if rep is not None),
              default=(None, None))
print('\nBEST REPAIRED OVERALL:', allbest)
allbestraw = max((raw, on) for on, (n, rf, raw, rep) in results.items())
print('BEST RAW OVERALL:', allbestraw)
print('DELIVERABLE BASELINE: 39026')
