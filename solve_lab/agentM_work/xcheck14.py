"""Widen the independent-verification scope of the enumeration.

T verified 9 of the 4,096 subsets of the 2^12 lattice against `checker.py`, all at the SHIPPED
granularity (nprobe=10). Round 14 prices at nprobe=80 and takes a max over greedy row orders,
which is a different code path through `tune`, so T's 9 points do not cover it.

This script materialises the full assignment my engine scores for a spread of subsets at the
p80 granularity and writes each to disk, so `checker.py` can be run on every one from outside
my parse. Subsets are chosen to SPAN the score range, not to cluster at the calibration point.
"""
import sys, os, json, time, itertools, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore                                            # noqa: E402

D4 = [642, 28730, 29854, 31864]
H16 = [642, 1844, 2892, 9629, 18253, 23642, 23754, 28355, 28730, 29305, 29854, 31864,
       34113, 35619, 37413, 37720]

# a spread: the witness, singletons, pairs, the |W|=3 best, supersets, a big support, and
# subsets DISJOINT from the witness (a region T did not sample at all).
CASES = [
    tuple(D4),                                        # the optimum
    (642,),                                           # singleton
    (28730,),                                         # singleton
    (23642,),                                         # singleton, non-witness
    (642, 28730),                                     # pair inside the witness
    (1844, 37413),                                    # pair disjoint from the witness
    (642, 28730, 29854),                              # witness minus one
    (642, 23642, 28730, 29854, 31864),                # |W|=5 superset
    (642, 2892, 23642, 28730, 29305, 29854, 31864),   # |W|=7 superset
    (1844, 2892, 9629, 18253, 23642),                 # |W|=5 disjoint from the witness
    tuple(sorted(H16[:10])),                          # |W|=10
    tuple(H16),                                       # the full support
]

out = []
for W in CASES:
    t0 = time.time()
    r = ieng.tune(list(W), nprobe=80, budget=180.0, want=True)
    freed, pin = ieng.site(list(W))
    if r.get('changes'):
        bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r['changes'], pin)
    else:
        bad, v = dict(ieng.BAD_UNC), list(ieng.V_UNC)
    sc = fscore.score(bad)
    fails = sorted(fscore.fails(bad))
    fn = f'xc14_{sc}_{len(W)}_{"-".join(map(str, W))}.json'
    if len(fn) > 110:
        fn = f'xc14_{sc}_{len(W)}_{abs(hash(W)) % 10**8}.json'
    json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0}, open(fn, 'w'))
    print(f'{fn}   engine score {sc}  ({len(fails)} failing)  {time.time()-t0:.1f}s', flush=True)
    out.append({'W': list(W), 'score': sc, 'nfail': len(fails), 'fails': fails[:12],
                'file': fn})
json.dump(out, open('xcheck14.json', 'w'), indent=1)
print('\nwrote xcheck14.json with', len(out), 'cases', flush=True)
