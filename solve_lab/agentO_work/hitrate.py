"""STEP 5 — THE RATE.

The admissible boundary changes delta form a coset delta0 + Lambda0 of a sublattice of Z^4
(the four configuration-tunable constants).  A configuration scan samples delta essentially at
random, so its hit rate is 1 / [Z^4 : Lambda0].  Measure the period in each direction with an
exact solvability oracle: the smallest t > 0 with delta0 + t*e_a still solvable.
"""
import sys, json, itertools, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/hitrate.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
SMALL = 2458959
BIGSMALL = 7376877


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const0, cols = G.build_model(R, Pv, G.V0)
Eqs, rows0 = G.eq_system(R, Pv, const0, cols)
b0 = {e: rows0[e][1] for e in Eqs}
Bcol = {a: {e: G.EQCO[a].get(e, 0) for e in Eqs if G.EQCO[a].get(e, 0)} for a in sorted(const0)}
TUN = [23616, 23618, 36660, 36662]


def solvable(delta):
    """delta: dict atom -> integer shift of that boundary constant."""
    rws, rhs = [], []
    for e in Eqs:
        r = {'z%d' % u: c for u, c in rows0[e][0].items()}
        rhs.append(b0[e] - sum(Bcol[a].get(e, 0) * d for a, d in delta.items()))
        rws.append(r)
    sol, msg, _ = sparse.solve_sparse(rws, rhs, names=list(range(len(Eqs))), verbose=False,
                                      maxcore=150, maxbits=10 ** 7, maxcorebits=10 ** 7)
    return sol is not None


# recover a particular admissible delta0 on the tunable support
rws, rhs = [], []
for e in Eqs:
    r = {'z%d' % u: c for u, c in rows0[e][0].items()}
    for a in TUN:
        c = Bcol[a].get(e, 0)
        if c:
            r['d%d' % a] = r.get('d%d' % a, 0) + c
    rws.append(r)
    rhs.append(b0[e])
sol, msg, _ = sparse.solve_sparse(rws, rhs, names=list(range(len(Eqs))), verbose=False,
                                  maxcore=150, maxbits=10 ** 7, maxcorebits=10 ** 7)
assert sol is not None, msg
# NOTE sign: the model above uses  A z + B d = b0, i.e. rhs shift of -B d
delta0 = {a: sol.get('d%d' % a, 0) for a in TUN}
say('a particular admissible boundary change delta0 (as a shift of const):')
for a in TUN:
    say('   const(a%-6d) shift = %d bits' % (a, abs(delta0[a]).bit_length()))
say('verify: solvable(delta0) = %s' % solvable(delta0))
say('verify: solvable(0)      = %s   (the witness itself)' % solvable({}))

say('\n--- period in each tunable direction (smallest t>0 with delta0 + t*e_a solvable)')
periods = {}
CANDS = [1, 2, 3, 9, 819653, SMALL, BIGSMALL, P, 3 * P, SMALL * P]
for a in TUN:
    per = None
    for t in CANDS:
        d = dict(delta0)
        d[a] = d[a] + t
        if solvable(d):
            per = t
            break
    periods[a] = per
    say('   direction const(a%-6d): period = %s'
        % (a, 'p' if per == P else ('3p' if per == 3 * P else
                                    ('2458959*p' if per == SMALL * P else str(per))) if per else '> all tested'))

say('\n--- index of the admissible lattice, assuming the measured periods are independent')
N = 1
unknown = False
for a in TUN:
    if periods[a] is None:
        unknown = True
    else:
        N *= periods[a]
if unknown:
    say('   at least one period exceeds every tested modulus; index is at least')
say('   [Z^4 : Lambda0] = %s' % (('%d bits ~ 2^%d' % (N.bit_length(), N.bit_length() - 1)) if N > 10 ** 12 else N))
say('   scan hit rate    = 1 / that = about 2^-%d' % (N.bit_length() - 1))
for M in (2800, 13884):
    say('   expected hits in %5d configurations = 2^-%d  (~10^-%d)'
        % (M, N.bit_length() - 1 - M.bit_length(),
           int((N.bit_length() - 1 - M.bit_length()) * 0.30103)))
json.dump({str(a): (str(periods[a]) if periods[a] else None) for a in TUN},
          open(OD + '/hitrate.json', 'w'))
say('DONE')
