"""Region growth around the 39,026 witness's residual.

A REGION is a set R of atoms that we allow to be nonzero.  Its PRIVATE variables are the
variables occurring in no atom outside R; those are free knobs that disturb nothing else.
The equations touched by R are E(R); every other equation stays satisfied because the private
variables cannot reach it.  Total failing equations = |E(R)| - maxsat(R), where maxsat(R) is
the largest subset of E(R) that some integer setting of the private variables satisfies.

The witness is R0 = {23616,23617,36659,36660,36661,36662,36663,36664}: 7 private variables,
|E| = 12, maxsat = 5, cost 7.  Growing R adds knobs and adds equations; we search for a region
with cost <= 6.
"""
import sys, json, itertools, collections, time, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
V0 = [0] * E.NV
for k, x in d.items():
    V0[int(k.split('_')[1])] = int(x)

R0 = [23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664]

# equations touching each atom
EQ_OF = collections.defaultdict(list)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        EQ_OF[a].append((e, c))
EQCO = {}
for a, lst in EQ_OF.items():
    EQCO[a] = dict(lst)


def private_vars(R):
    Rs = set(R)
    P = set()
    for a in R:
        for u in H.avars[a]:
            if set(H.occ[u]) <= Rs:
                P.add(u)
    return sorted(P)


def atom_val(w, a):
    ns = {'v': w, '__builtins__': {}}
    return eval(H.acodes[a], ns)


def build_model(R, P, base):
    """Residual of each atom in R as an affine function of the private vars P.
       Returns (const {a:val}, cols {u:{a:coef}}) or None if any dependence is nonlinear."""
    w = list(base)
    const = {a: atom_val(w, a) for a in R}
    cols = {}
    for u in P:
        o = w[u]
        w[u] = o + 1
        d1 = {a: atom_val(w, a) - const[a] for a in R}
        w[u] = o + 2
        d2 = {a: atom_val(w, a) - const[a] for a in R}
        w[u] = o
        for a in R:
            if d2[a] != 2 * d1[a]:
                return None
        cols[u] = {a: c for a, c in d1.items() if c}
    return const, cols


def eq_system(R, P, const, cols):
    Eqs = sorted({e for a in R for e in EQCO[a]})
    rows = {}
    for e in Eqs:
        r = {}
        for u in P:
            s = 0
            for a, c in cols[u].items():
                s += EQCO[a].get(e, 0) * c
            if s:
                r[u] = s
        rhs = -sum(EQCO[a].get(e, 0) * const[a] for a in R)
        rows[e] = (r, rhs)
    return Eqs, rows


def solvable(S, rows):
    sol, msg, _ = sparse.solve_sparse([rows[e][0] for e in S], [rows[e][1] for e in S],
                                      names=list(S), verbose=False, maxcore=80,
                                      maxbits=10 ** 7, maxcorebits=10 ** 7)
    return sol


def maxsat(Eqs, rows, exhaustive_upto=15, tries=400, seed=3):
    n = len(Eqs)
    if n <= exhaustive_upto:
        for k in range(n, 0, -1):
            for S in itertools.combinations(Eqs, k):
                sol = solvable(S, rows)
                if sol is not None:
                    return k, list(S), sol
        return 0, [], {}
    rnd = random.Random(seed)
    best = (0, [], {})
    for t in range(tries):
        order = list(Eqs)
        if t:
            rnd.shuffle(order)
        keep = []
        sol = {}
        for e in order:
            cand = keep + [e]
            s = solvable(cand, rows)
            if s is not None:
                keep = cand
                sol = s
        if len(keep) > best[0]:
            best = (len(keep), sorted(keep), sol)
    return best


def evaluate(R, base=V0, verbose=True):
    R = sorted(set(R))
    P = private_vars(R)
    m = build_model(R, P, base)
    if m is None:
        return None
    const, cols = m
    Eqs, rows = eq_system(R, P, const, cols)
    k, S, sol = maxsat(Eqs, rows)
    cost = len(Eqs) - k
    if verbose:
        print(f'  R({len(R)} atoms) P={len(P)} |E|={len(Eqs)} maxsat={k} COST={cost}', flush=True)
    return cost, len(Eqs), k, S, sol, P, R


def realise(R, P, sol, base=V0):
    w = list(base)
    for u in P:
        w[u] = sol.get(u, 0)
    return w


if __name__ == '__main__':
    print('baseline region:')
    r = evaluate(R0)
    cost0 = r[0]
    print('  cost', cost0, '(expect 7)')
    w = realise(r[6], r[5], r[4])
    print('  exact check:', len(E.eqfails(E.badatoms(w))), 'failing equations')

    # candidate atoms to absorb: atoms sharing a variable with R0
    cand = set()
    for a in R0:
        for u in H.avars[a]:
            cand |= set(H.occ[u])
    cand -= set(R0)
    print(f'\n{len(cand)} adjacent atoms; single-atom growth:', flush=True)
    results = []
    t0 = time.time()
    for a in sorted(cand):
        R = R0 + [a]
        try:
            r = evaluate(R, verbose=False)
        except Exception as ex:
            print(f'  +a{a}: ERR {type(ex).__name__}', flush=True)
            continue
        if r is None:
            print(f'  +a{a}: nonlinear in a private var', flush=True)
            continue
        cost, ne, k, S, sol, P, RR = r
        results.append((cost, a, ne, k, len(P)))
        flag = ' ***BETTER***' if cost < cost0 else ''
        print(f'  +a{a}: P={len(P)} |E|={ne} maxsat={k} COST={cost}{flag}', flush=True)
        if cost < cost0:
            w = realise(RR, P, sol)
            nf = len(E.eqfails(E.badatoms(w)))
            print(f'      exact re-evaluation: {nf} failing -> score {39033-nf}', flush=True)
            if nf < 7:
                json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0},
                          open(f'{OD}/grow_{a}_{39033-nf}.json', 'w'))
                print('      *** WROTE improvement', flush=True)
    results.sort()
    print(f'\nbest single-atom growths ({time.time()-t0:.0f}s):', results[:12])
    json.dump([[int(x) for x in r] for r in results], open(f'{OD}/grow1.json', 'w'))
