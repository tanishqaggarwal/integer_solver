"""Given a handle set X whose footprint is a small equation set E, maximise the
number of equations in E satisfied (allowing currently-satisfied ones to break).
Exact over Z."""
import sys, itertools, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, push, linmodel, intsolve, modp


def footprint(X):
    E = set()
    for x in X:
        for a in lib.var_atoms[x]:
            E.update(lib.atom2eq.get(a, ()))
    return sorted(E)


def optimise(V, X, verbose=True, kmin=1):
    assert not push.conflicts(X), 'nonlinear conflict'
    rows, rhs, issq, eqs = linmodel.build(X, V)
    E = [e for e in eqs if any(rows[e])]
    n = len(X)
    cur_fail = [e for e in E if rhs[e] != 0]
    if verbose:
        print(f'|X|={n} |E|={len(E)} currently failing in E: {len(cur_fail)} {cur_fail}')
    for k in range(len(E), kmin - 1, -1):
        for sub in itertools.combinations(E, k):
            M = [rows[e] for e in sub]
            b = [-rhs[e] for e in sub]
            if not modp.solvable(M, b, n):
                continue
            D, info = intsolve.solve_int(M, b, n)
            if D is not None:
                return k, sub, D, E
        if verbose:
            print(f'  k={k}: none')
    return 0, (), None, E


if __name__ == '__main__':
    STATE = sys.argv[1]
    X = [int(z) for z in sys.argv[2].split(',')]
    V = lib.H.load_assignment(STATE)
    F0 = lib.true_fails(V)
    print('state', STATE, 'failing', len(F0), F0, flush=True)
    t0 = time.time()
    k, sub, D, E = optimise(V, X)
    print(f'max satisfiable in E: {k}/{len(E)}  ({time.time()-t0:.0f}s)')
    if D:
        W = list(V)
        for i, x in enumerate(X):
            W[x] += D[i]
        f = lib.true_fails(W)
        print('VERIFIED failing:', len(f), f)
        if len(f) < len(F0):
            out = sys.argv[3] if len(sys.argv) > 3 else f'sa/opt_{len(f)}.json'
            lib.H.save_assignment(W, out)
            with open('sa/hits.jsonl', 'a') as fh:
                fh.write(json.dumps({'method': 'opt', 'state_in': STATE, 'X': X,
                                     'D': [str(z) for z in D], 'kept': list(sub),
                                     'nfail': len(f), 'failing': f, 'out': out}) + '\n')
            print('saved', out)
