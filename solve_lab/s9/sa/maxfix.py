"""Maximise the number of simultaneously-zeroable failing equations from a given
state, over a chosen perturbation variable set, exactly over Z."""
import sys, itertools, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve

STATE = sys.argv[1] if len(sys.argv) > 1 else '../best/new_instance_partial_39022.json'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'sa/maxfix_out.json'

V = lib.H.load_assignment(STATE)
F = lib.true_fails(V)
Fs = set(F)
print('state', STATE, 'failing', len(F), F, flush=True)

order = [9413, 17325, 642, 28730, 29854, 31864, 23754, 1844, 35619, 29305, 9629,
         1329, 10903, 6947, 21574, 33168, 1613, 950, 6090, 2892, 15120, 28355,
         18253, 35531, 10422, 37720, 23642, 11099, 37413]
X = []
for x in order:
    if push.conflicts(X + [x]):
        continue
    X.append(x)
print('|X| =', len(X), X, flush=True)

rows, rhs, issq, eqs = linmodel.build(X, V)
touched = [e for e in eqs if any(rows[e])]
must = [e for e in touched if e not in Fs]
cand = [e for e in touched if e in Fs]
print('touched', len(touched), 'must', len(must), 'cand', len(cand), cand, flush=True)

n = len(X)
t0 = time.time()
best = None
for k in range(len(cand), 0, -1):
    tried = 0
    for sub in itertools.combinations(cand, k):
        M = [rows[e] for e in must] + [rows[e] for e in sub]
        b = [0] * len(must) + [-rhs[e] for e in sub]
        D, info = intsolve.solve_int(M, b, n)
        tried += 1
        if D is not None:
            best = (k, sub, D)
            break
    print(f'k={k}: {tried} subsets, {"HIT" if best else "none"}  {time.time()-t0:.0f}s', flush=True)
    if best:
        break

if best:
    k, sub, D = best
    W = list(V)
    for i, x in enumerate(X):
        W[x] += D[i]
    f = lib.true_fails(W)
    print('fixed subset', sub)
    print('VERIFIED failing:', len(f), f, flush=True)
    if len(f) < len(F):
        lib.H.save_assignment(W, OUT)
        print('saved', OUT)
        with open('sa/hits.jsonl', 'a') as fh:
            fh.write(json.dumps({'method': 'maxfix', 'state_in': STATE,
                                 'X': X, 'D': [str(d) for d in D],
                                 'fixed': list(sub), 'nfail': len(f),
                                 'failing': f, 'out': OUT}) + '\n')
else:
    print('no improvement possible with this variable set')
