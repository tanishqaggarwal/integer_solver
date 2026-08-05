"""Stochastic + beam search over perturbation-variable subsets.
usage: python3 sa/search2.py <state.json> <mode: beam|rand> <budget_s> <prefix> [seed]
"""
import sys, itertools, pickle, time, json, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve, fastmax

_d = pickle.load(open('atoms.pkl', 'rb'))
ET = _d['eq_terms']

STATE = sys.argv[1]
MODE = sys.argv[2]
BUDGET = float(sys.argv[3])
PREF = sys.argv[4]
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 0
random.seed(SEED)

V = lib.H.load_assignment(STATE)
F = lib.true_fails(V)
Fs = set(F)
print('state', STATE, 'failing', len(F), F, flush=True)

atoms = set()
for e in F:
    for c, a in ET[e][2]:
        atoms.add(a)
cv = {}
for a in sorted(atoms):
    for m in push.apoly(a):
        for x in m:
            s = set()
            for at in lib.var_atoms[x]:
                s.update(lib.atom2eq.get(at, ()))
            cv[x] = min(cv.get(x, 10**9), len(s - Fs))
CAND = [x for x, k in sorted(cv.items(), key=lambda t: (t[1], t[0]))]
print('candidate vars:', len(CAND), flush=True)

CACHE = {}
NEVAL = 0


def evaluate(X):
    global NEVAL
    key = tuple(X)
    if key in CACHE:
        return CACHE[key]
    if push.conflicts(X):
        CACHE[key] = None
        return None
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [rows[e] + [0] for e in touched if e not in Fs]
    cand = [(e, rows[e] + [-rhs[e]]) for e in touched if e in Fs]
    NEVAL += 1
    if not cand:
        r = (0, (), len(must))
    else:
        k, sub = fastmax.max_subset(must, cand, len(X))
        r = (k, sub, len(must))
    CACHE[key] = r
    return r


def verify(X, sub):
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [e for e in touched if e not in Fs]
    n = len(X)
    M = [rows[e] for e in must] + [rows[e] for e in sub]
    b = [0] * len(must) + [-rhs[e] for e in sub]
    D, info = intsolve.solve_int(M, b, n)
    if D is None:
        return None, info
    W = list(V)
    for i, x in enumerate(X):
        W[x] += D[i]
    return W, D


def record(X, sub, W, D):
    f = lib.true_fails(W)
    print('  EXACT-Z verified: failing', len(f), f, flush=True)
    if len(f) < len(F):
        out = f'{PREF}_{len(f)}.json'
        lib.H.save_assignment(W, out)
        with open('sa/hits.jsonl', 'a') as fh:
            fh.write(json.dumps({'method': MODE, 'state_in': STATE, 'X': X,
                                 'D': [str(z) for z in D], 'fixed': list(sub),
                                 'nfail': len(f), 'failing': f, 'out': out}) + '\n')
        print('  saved', out, flush=True)
    return len(f)


t0 = time.time()
best = (0, None, None)

if MODE == 'beam':
    BW = 12
    beam = [[]]
    d = 0
    while time.time() - t0 < BUDGET:
        nxt = []
        for X in beam:
            for x in CAND:
                if x in X:
                    continue
                Y = sorted(X + [x])
                r = evaluate(Y)
                if r is None:
                    continue
                nxt.append((r[0], r[2], Y, r[1]))
            if time.time() - t0 > BUDGET:
                break
        if not nxt:
            break
        nxt.sort(key=lambda t: (-t[0], t[1], len(t[2])))
        d += 1
        print(f'depth {d}: {len(nxt)} states  best k={nxt[0][0]} must={nxt[0][1]} '
              f'X={nxt[0][2]} sub={nxt[0][3]}  evals={NEVAL} {time.time()-t0:.0f}s', flush=True)
        if nxt[0][0] > best[0]:
            best = (nxt[0][0], nxt[0][2], nxt[0][3])
            W, D = verify(best[1], best[2])
            if W is not None:
                record(best[1], best[2], W, D)
            else:
                print('  (mod-P only; no integer solution:', str(D)[:80], ')', flush=True)
        beam = [t[2] for t in nxt[:BW]]
else:
    while time.time() - t0 < BUDGET:
        sz = random.randint(2, 14)
        X = sorted(random.sample(CAND, sz))
        r = evaluate(X)
        if r is None or r[0] <= best[0]:
            continue
        best = (r[0], X, r[1])
        print(f'rand: k={r[0]} must={r[2]} X={X} sub={r[1]} evals={NEVAL} '
              f'{time.time()-t0:.0f}s', flush=True)
        W, D = verify(X, r[1])
        if W is not None:
            record(X, r[1], W, D)
        else:
            print('  (mod-P only; no integer solution)', flush=True)

print('done. evals =', NEVAL, 'best k =', best[0], flush=True)
print('best X =', best[1], 'sub =', best[2], flush=True)
