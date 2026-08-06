"""Fast randomized greedy for the maximum feasible subsystem (mod P), using
streaming Gaussian elimination, then exact Z verification.

usage: python3 sa/greedy2.py <state.json> <nvars> <seconds> <prefix> <seed>
"""
import sys, random, time, json, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve, bighandle

P = lib.P
STATE = sys.argv[1]
NVARS = int(sys.argv[2])
BUDGET = float(sys.argv[3])
PREF = sys.argv[4]
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 0
random.seed(SEED)

V = lib.H.load_assignment(STATE)
F0 = lib.true_fails(V)
Fs = set(F0)
print('state', STATE, 'failing', len(F0), F0, flush=True)

X = bighandle.X if False else None
# rebuild the handle set here (bighandle does it at import time for its own argv)
import pickle as _p
ET = _p.load(open('atoms.pkl', 'rb'))['eq_terms']


def eff_foot(x):
    return bighandle.eff_foot(x, V)


L1 = bighandle.vars_of_eqs(F0)
seedE = set(F0)
for x in L1:
    seedE |= eff_foot(x)
L2 = bighandle.vars_of_eqs(seedE)
CANDS = sorted(L1 | L2)
scored = []
for x in CANDS:
    E = eff_foot(x)
    if not E:
        continue
    scored.append((len(E - Fs), -len(E & Fs), x))
scored.sort()
X = []
for _, _, x in scored:
    if len(X) >= NVARS:
        break
    if bighandle.ok_to_add(X, x):
        X.append(x)
print('|X| =', len(X), flush=True)

rows, rhs, issq, eqs = linmodel.build(X, V)
E = [e for e in eqs if any(rows[e])]
n = len(X)
R = {e: [c % P for c in rows[e]] + [(-rhs[e]) % P] for e in E}
print('|E| =', len(E), 'failing inside E:', len([e for e in E if rhs[e] != 0]), flush=True)


def greedy(order):
    piv = []          # list of (col, row) already reduced
    S = []
    for e in order:
        r = R[e][:]
        for c, pr in piv:
            if r[c]:
                f = r[c]
                r = [(a - f * b) % P for a, b in zip(r, pr)]
        c0 = None
        for c in range(n):
            if r[c]:
                c0 = c; break
        if c0 is None:
            if r[n] == 0:
                S.append(e)          # redundant but consistent
            continue                  # inconsistent -> drop e
        inv = pow(r[c0], P - 2, P)
        r = [(z * inv) % P for z in r]
        piv.append((c0, r))
        S.append(e)
    return S


best = (0, None)
t0 = time.time()
it = 0
while time.time() - t0 < BUDGET:
    it += 1
    order = E[:]
    random.shuffle(order)
    S = greedy(order)
    if len(S) > best[0]:
        best = (len(S), sorted(S))
        print(f'iter {it}: |S|={len(S)} unsat={len(E)-len(S)} {time.time()-t0:.0f}s', flush=True)
print(f'{it} passes; best |S|={best[0]}/{len(E)} -> {len(E)-best[0]} unsatisfied', flush=True)

S = best[1]
if S:
    D, info = intsolve.solve_int([rows[e] for e in S], [-rhs[e] for e in S], n)
    tries = 0
    while D is None and len(S) > 1 and tries < 40:
        S = S[:-1]
        tries += 1
        D, info = intsolve.solve_int([rows[e] for e in S], [-rhs[e] for e in S], n)
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
                fh.write(json.dumps({'method': 'greedy2', 'state_in': STATE, 'X': X,
                                     'D': [str(z) for z in D], 'nfail': len(f),
                                     'failing': f, 'out': out}) + '\n')
            print('saved', out, flush=True)
    else:
        print('no integer solution found', flush=True)
