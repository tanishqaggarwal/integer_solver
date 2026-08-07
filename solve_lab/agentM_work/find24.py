"""Locate the 39,024 points the cofactor knob set produced, and materialise them.

`enumcof.py` at |W| = 5 reports `39024:3` -- a score that appears NOWHERE in any enumeration
run with handle-only knobs (the no-cofactor spectrum at |W| = 5 goes 39,026 then straight to
39,022).  39,024 is the best score ever reached at a placement OTHER than the witness, beating
the previous non-witness record of 39,023.  The checkpoint stores distributions only, so this
re-scans the sizes named on the command line and records the actual subsets.
"""
import sys, os, json, time, itertools, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore                                            # noqa: E402
import enumcof_lib as L                                        # noqa: E402

SIZES = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '4,5').split(',')]
THRESH = int(sys.argv[2]) if len(sys.argv) > 2 else 39023

HL = L.HL
hits = []
t0 = time.time()
n = 0
for k in SIZES:
    for W in itertools.combinations(HL, k):
        n += 1
        try:
            r = L.tune_ext(list(W))
        except Exception:
            continue
        if r and r['score'] >= THRESH:
            hits.append((r['score'], W))
        if n % 1000 == 0:
            print(f'  [{n}] {time.time()-t0:.0f}s  hits {len(hits)}', flush=True)

hits.sort(reverse=True)
print(f'\n{len(hits)} subsets at score >= {THRESH} over |W| in {SIZES}  ({time.time()-t0:.0f}s)',
      flush=True)
out = []
for sc, W in hits[:20]:
    r = L.tune_ext(list(W))
    if r['changes']:
        bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r['changes'], r['pin'])
    else:
        bad, v = dict(ieng.BAD_UNC), list(ieng.V_UNC)
    s2 = fscore.score(bad)
    fn = f'M_cof24_{s2}_{"-".join(map(str, W))}.json'
    if len(fn) > 110:
        fn = f'M_cof24_{s2}_{abs(hash(W)) % 10**8}.json'
    json.dump({f"x_{i}": int(v[i]) for i in range(ieng.NV) if v[i] != 0}, open(fn, 'w'))
    print(f'  {s2}  W={W}  -> {fn}', flush=True)
    out.append({'score': s2, 'W': list(W), 'file': fn,
                'fails': sorted(fscore.fails(bad))})
json.dump(out, open('find24.json', 'w'), indent=1)
