"""Beam search over perturbation-variable subsets, maximising the number of
failing equations that can be zeroed (mod-P screen, then exact Z verification).

usage: python3 sa/beam.py <state.json> <beamwidth> <depth> [outprefix]
"""
import sys, itertools, pickle, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve, modp

_d = pickle.load(open('atoms.pkl', 'rb'))
ET = _d['eq_terms']

STATE = sys.argv[1]
BW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
DEPTH = int(sys.argv[3]) if len(sys.argv) > 3 else 16
PREF = sys.argv[4] if len(sys.argv) > 4 else 'sa/beam'

V = lib.H.load_assignment(STATE)
F = lib.true_fails(V)
Fs = set(F)
print('state', STATE, 'failing', len(F), F, flush=True)

# candidate variables: those appearing in atoms of the failing equations,
# plus one level out (atoms sharing an equation with them) -- ranked by footprint.
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


def evaluate(X):
    """(best k solvable mod P, must size, rows/rhs/must/cand)."""
    if push.conflicts(X):
        return None
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [e for e in touched if e not in Fs]
    cand = [e for e in touched if e in Fs]
    n = len(X)
    for k in range(len(cand), 0, -1):
        for sub in itertools.combinations(cand, k):
            M = [rows[e] for e in must] + [rows[e] for e in sub]
            b = [0] * len(must) + [-rhs[e] for e in sub]
            if modp.solvable(M, b, n):
                return k, sub, len(must), (rows, rhs, must, cand)
    return 0, (), len(must), (rows, rhs, must, cand)


beam = [([], 0, 0)]
seen = set()
best_overall = (0, None, None)
t0 = time.time()
for d in range(DEPTH):
    nxt = []
    for X, k0, m0 in beam:
        for x in CAND:
            if x in X:
                continue
            Y = sorted(X + [x])
            key = tuple(Y)
            if key in seen:
                continue
            seen.add(key)
            r = evaluate(Y)
            if r is None:
                continue
            k, sub, nmust, data = r
            nxt.append((Y, k, nmust, sub))
    if not nxt:
        break
    nxt.sort(key=lambda t: (-t[1], t[2], len(t[0])))
    print(f'depth {d+1}: {len(nxt)} states, best k={nxt[0][1]} must={nxt[0][2]} X={nxt[0][0]} '
          f'sub={nxt[0][3]}  {time.time()-t0:.0f}s', flush=True)
    if nxt[0][1] > best_overall[0]:
        best_overall = (nxt[0][1], nxt[0][0], nxt[0][3])
    beam = [(t[0], t[1], t[2]) for t in nxt[:BW]]

print('BEST mod-P:', best_overall, flush=True)

# exact Z verification on the top candidates
k, X, sub = best_overall
if X:
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [e for e in touched if e not in Fs]
    cand = [e for e in touched if e in Fs]
    n = len(X)
    for kk in range(k, 0, -1):
        done = False
        for s in itertools.combinations(cand, kk):
            M = [rows[e] for e in must] + [rows[e] for e in s]
            b = [0] * len(must) + [-rhs[e] for e in s]
            if not modp.solvable(M, b, n):
                continue
            D, info = intsolve.solve_int(M, b, n)
            if D is None:
                continue
            W = list(V)
            for i, xx in enumerate(X):
                W[xx] += D[i]
            f = lib.true_fails(W)
            print('EXACT: fixing', s, '-> failing', len(f), f, flush=True)
            if len(f) < len(F):
                out = f'{PREF}_{len(f)}.json'
                lib.H.save_assignment(W, out)
                with open('sa/hits.jsonl', 'a') as fh:
                    fh.write(json.dumps({'method': 'beam', 'state_in': STATE, 'X': X,
                                         'D': [str(z) for z in D], 'fixed': list(s),
                                         'nfail': len(f), 'failing': f, 'out': out}) + '\n')
                print('saved', out)
                done = True
            break
        if done:
            break
