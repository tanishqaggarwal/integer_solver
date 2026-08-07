"""Exact optimisation over the SEVEN region-private integers.

The 39,026 witness's whole residual lives in 8 atoms touched by exactly 12 equations, and
seven variables occur in NO atom outside that region:
    x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864
Every one of the 8 residuals is an affine function of those seven, so the region reduces to
12 linear equations in 7 integer unknowns.  Maximise the number satisfied; the witness
satisfies 5 of 12 (7 failing).  Anything better than 5 beats 39,026.
"""
import sys, json, itertools, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v = [0] * E.NV
for k, x in d.items():
    v[int(k.split('_')[1])] = int(x)

ATOMS = [23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664]
bad = E.badatoms(v)
touch = collections.defaultdict(dict)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        if a in ATOMS:
            touch[e][a] = c
EQ = sorted(touch)
assert set(bad) == set(ATOMS), sorted(bad)

K1 = v[7068] - v[2099]
K2 = 5113045 * (v[7075] * v[9118])
J = v[7075] * v[8731]
x17499, x22665, x28961, x28599, x28730 = v[17499], v[22665], v[28961], v[28599], v[28730]

U = [642, 9413, 1329, 29854, 10903, 31864, 17325]


def row_of(e):
    c = lambda a: touch[e].get(a, 0)
    r = {
        642:   -7376877 * c(23616) + c(36664),
        9413:  -x17499 * c(23617),
        1329:  -x22665 * c(36659),
        29854: c(36659) - c(36660),
        10903: -x28961 * c(36661),
        31864: c(36661) + c(36663),
        17325: -x28599 * c(36664),
    }
    const = K1 * c(23616) + x28730 * c(23617) + K2 * c(36660) + J * c(36662)
    return {k: x for k, x in r.items() if x}, -const     # row . z = rhs


ROWS = {e: row_of(e) for e in EQ}

# sanity: the witness's own values must satisfy exactly the 5 currently-satisfied equations
z0 = {u: v[u] for u in U}
sat0 = [e for e in EQ if sum(cc * z0[k] for k, cc in ROWS[e][0].items()) == ROWS[e][1]]
cur = [e for e in EQ if sum(cc * bad.get(a, 0) for a, cc in touch[e].items()) == 0]
print('witness satisfies (via reduced model):', sat0)
print('witness satisfies (direct):           ', cur)
assert sat0 == cur, 'reduced model does not reproduce the witness'
print('reduced model VERIFIED against the witness\n')

best = None
for k in range(len(EQ), 0, -1):
    if best is not None:
        break
    found = []
    for S in itertools.combinations(EQ, k):
        rows = [ROWS[e][0] for e in S]
        rhs = [ROWS[e][1] for e in S]
        sol, msg, _ = sparse.solve_sparse(rows, rhs, names=list(S), verbose=False,
                                          maxcore=64, maxbits=10 ** 7, maxcorebits=10 ** 7)
        if sol is None:
            continue
        z = {u: sol.get(u, 0) for u in U}
        w = list(v)
        for u in U:
            w[u] = z[u]
        ba = E.badatoms(w)
        ff = E.eqfails(ba)
        found.append((len(ff), S, z, sorted(ba)))
    if found:
        found.sort(key=lambda t: t[0])
        print(f'size {k}: {len(found)} integrally solvable subsets, best fails={found[0][0]}')
        best = (k, found)
    else:
        print(f'size {k}: none integrally solvable')

k, found = best
nf, S, z, ba = found[0]
print('\nBEST: satisfied subset', S)
print('  fails =', nf, ' score =', 39033 - nf, ' bad atoms =', ba)
w = list(v)
for u in U:
    w[u] = z[u]
out = f'{OD}/region_opt_{39033-nf}.json'
json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0}, open(out, 'w'))
print('wrote', out)
for u in U:
    print(f'  x_{u}: {abs(v[u]).bit_length()}b -> {abs(z[u]).bit_length()}b')
