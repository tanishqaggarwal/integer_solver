"""STEP 3 — invert rather than sample.

Integral solvability of the region system A z = b is a coset condition on b, and b is linear
in the region's boundary constants const_a (a in R) -- exactly what a configuration change
moves.  So instead of sampling configurations, put the boundary perturbation delta into the
unknown vector and solve

        A z  +  B delta  =  b0        over Z,

where B is the map (boundary perturbation) -> (rhs).  A solution says which boundary constants
must change and by how much; no solution with a given support says that support can never work,
however many configurations are sampled.
"""
import sys, json, itertools, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/invert.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const0, cols = G.build_model(R, Pv, G.V0)
Eqs, rows0 = G.eq_system(R, Pv, const0, cols)
say('region R = %s' % sorted(R))
say('private unknowns z: %s' % Pv)
say('equations: %s' % Eqs)

# rhs as a function of boundary perturbation: b_e(delta) = b0_e - sum_a c_{e,a} delta_a
Bcol = {}
for a in sorted(const0):
    Bcol[a] = {e: G.EQCO[a].get(e, 0) for e in Eqs if G.EQCO[a].get(e, 0)}
b0 = {e: rows0[e][1] for e in Eqs}


def try_support(S, tag):
    """Solve A z + sum_{a in S} (column of a) * delta_a = b0 over Z."""
    rws, rhs = [], []
    for e in Eqs:
        r = {'z%d' % u: c for u, c in rows0[e][0].items()}
        for a in S:
            c = Bcol[a].get(e, 0)
            if c:
                r['d%d' % a] = r.get('d%d' % a, 0) + c
        rws.append(r)
        rhs.append(b0[e])
    sol, msg, _ = sparse.solve_sparse(rws, rhs, names=list(range(len(Eqs))), verbose=False,
                                      maxcore=150, maxbits=10 ** 7, maxcorebits=10 ** 7)
    if sol is None:
        say('  %-34s UNSOLVABLE over Z: %s' % (tag, msg[:70]))
        return None
    d = {a: sol.get('d%d' % a, 0) for a in S}
    say('  %-34s SOLVABLE.  required boundary change:' % tag)
    for a in S:
        v = d[a]
        say('        const(a%-6d) += %s  (%d bits)%s'
            % (a, (str(v)[:26] + '..') if abs(v) > 10 ** 20 else v, abs(v).bit_length(),
               '   [== 0, no change needed]' if v == 0 else ''))
    return d


say('\n--- baseline: no boundary change allowed')
try_support([], 'delta = 0')

say('\n--- single boundary constant free')
ok1 = []
for a in sorted(const0):
    r = try_support([a], 'only const(a%d) free' % a)
    if r is not None:
        ok1.append(a)
say('  single-constant supports that work: %s' % ok1)

say('\n--- all boundary constants free')
allfree = try_support(sorted(const0), 'all 9 constants free')

if not ok1 and allfree is None:
    say('\n  => no change to the region boundary can make this region integrally solvable.')
elif ok1:
    say('\n  => the obstruction is invertible from a SINGLE boundary constant.')

say('\n--- minimal supports (pairs), if singles failed')
if not ok1:
    good = []
    for S in itertools.combinations(sorted(const0), 2):
        r = try_support(list(S), 'const(a%d,a%d) free' % S)
        if r is not None:
            good.append(S)
    say('  working pairs: %s' % good)
say('DONE')
