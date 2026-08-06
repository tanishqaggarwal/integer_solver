"""Scalable version of sa/opt.py: randomized greedy maximisation of the number of
satisfied equations in the footprint of a handle set (mod-P feasibility, then
exact Z verification).

usage: python3 sa/greedyopt.py <state.json> <vars csv|auto> <seconds> <outprefix>
"""
import sys, random, time, json, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve, modp

P = lib.P
STATE = sys.argv[1]
VARSPEC = sys.argv[2]
BUDGET = float(sys.argv[3])
PREF = sys.argv[4]
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 0
random.seed(SEED)

V = lib.H.load_assignment(STATE)
F0 = lib.true_fails(V)
print('state', STATE, 'failing', len(F0), F0, flush=True)

if VARSPEC == 'auto':
    import pickle
    ET = pickle.load(open('atoms.pkl', 'rb'))['eq_terms']
    Fs = set(F0)
    atoms = set()
    for e in F0:
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
    ranked = [x for x, k in sorted(cv.items(), key=lambda t: (t[1], t[0]))]
    X = []
    for x in ranked:
        if push.conflicts(X + [x]):
            continue
        X.append(x)
else:
    X = [int(z) for z in VARSPEC.split(',')]
    assert not push.conflicts(X), 'nonlinear conflict'

print('|X| =', len(X), X, flush=True)
rows, rhs, issq, eqs = linmodel.build(X, V)
E = [e for e in eqs if any(rows[e])]
n = len(X)
print('|E| =', len(E), 'currently failing inside E:',
      len([e for e in E if rhs[e] != 0]), flush=True)


def feasible(S):
    M = [rows[e] for e in S]
    b = [-rhs[e] for e in S]
    return modp.solvable(M, b, n)


best = (0, None)
t0 = time.time()
it = 0
while time.time() - t0 < BUDGET:
    it += 1
    order = E[:]
    random.shuffle(order)
    S = []
    for e in order:
        if feasible(S + [e]):
            S.append(e)
    if len(S) > best[0]:
        best = (len(S), sorted(S))
        print(f'iter {it}: |S|={len(S)}  ({len(E)-len(S)} unsatisfied in E)  '
              f'{time.time()-t0:.0f}s', flush=True)
print(f'{it} greedy passes; best mod-P |S| = {best[0]} of {len(E)}', flush=True)

S = best[1]
if S:
    M = [rows[e] for e in S]
    b = [-rhs[e] for e in S]
    D, info = intsolve.solve_int(M, b, n)
    if D is None:
        print('no integer solution for the best mod-P set:', str(info)[:100], flush=True)
        # fall back: shrink until integral
        cur = S[:]
        while cur and D is None:
            cur = cur[:-1]
            D, info = intsolve.solve_int([rows[e] for e in cur], [-rhs[e] for e in cur], n)
        S = cur
    if D is not None:
        W = list(V)
        for i, x in enumerate(X):
            W[x] += D[i]
        f = lib.true_fails(W)
        print('EXACT-Z verified failing:', len(f), f, flush=True)
        if len(f) < len(F0):
            out = f'{PREF}_{len(f)}.json'
            lib.H.save_assignment(W, out)
            with open('sa/hits.jsonl', 'a') as fh:
                fh.write(json.dumps({'method': 'greedyopt', 'state_in': STATE, 'X': X,
                                     'D': [str(z) for z in D], 'kept': list(S),
                                     'nfail': len(f), 'failing': f, 'out': out}) + '\n')
            print('saved', out, flush=True)
