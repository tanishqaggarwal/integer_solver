"""Greedy max-feasible-subsystem with a SEEDED handle set (e.g. the p-wire var
x_17499 plus its trackers), then filled with the ranked local handles.

usage: python3 sa/wire.py <state.json> <seedvars csv> <nvars> <seconds> <prefix> <seed>
"""
import sys, random, time, json, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, linmodel, intsolve, bighandle

P = lib.P
STATE = sys.argv[1]
SEEDV = [int(z) for z in sys.argv[2].split(',')] if sys.argv[2] else []
NVARS = int(sys.argv[3])
BUDGET = float(sys.argv[4])
PREF = sys.argv[5]
random.seed(int(sys.argv[6]) if len(sys.argv) > 6 else 0)

V = lib.H.load_assignment(STATE)
F0 = lib.true_fails(V)
Fs = set(F0)
print('state', STATE, 'failing', len(F0), F0, flush=True)

X = []
for x in SEEDV:
    if bighandle.ok_to_add(X, x):
        X.append(x)
    else:
        print('  seed x_%d rejected (nonlinear conflict)' % x, flush=True)
print('seeded |X| =', len(X), X, flush=True)

L1 = bighandle.vars_of_eqs(F0)
seedE = set(F0)
for x in L1:
    seedE |= bighandle.eff_foot(x, V)
CANDS = sorted(L1 | bighandle.vars_of_eqs(seedE))
scored = []
for x in CANDS:
    if x in X:
        continue
    E = bighandle.eff_foot(x, V)
    if not E:
        continue
    scored.append((len(E - Fs), -len(E & Fs), x))
scored.sort()
for _, _, x in scored:
    if len(X) >= NVARS:
        break
    if bighandle.ok_to_add(X, x):
        X.append(x)
print('|X| =', len(X), flush=True)

t0 = time.time()
rows, rhs, issq, eqs = linmodel.build(X, V)
E = [e for e in eqs if any(rows[e])]
n = len(X)
R = {e: [c % P for c in rows[e]] + [(-rhs[e]) % P] for e in E}
print('|E| =', len(E), 'failing inside E:', len([e for e in E if rhs[e] != 0]),
      f'(build {time.time()-t0:.0f}s)', flush=True)


def greedy(order):
    piv = []
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
                S.append(e)
            continue
        inv = pow(r[c0], P - 2, P)
        piv.append((c0, [(z * inv) % P for z in r]))
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
        S = S[:-1]; tries += 1
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
                fh.write(json.dumps({'method': 'wire', 'state_in': STATE, 'X': X,
                                     'D': [str(z) for z in D], 'nfail': len(f),
                                     'failing': f, 'out': out}) + '\n')
            print('saved', out, flush=True)
    else:
        print('no integer solution', flush=True)
