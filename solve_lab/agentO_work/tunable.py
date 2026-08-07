"""STEP 4 — which boundary constants can a CONFIGURATION actually move?

Region atoms and their external ingredients (everything else in them is private):
  a23616 = x_7068 - x_2099 - 7376877*x_642      external: K1 = x_7068 - x_2099      TUNABLE
  a23617 = x_28730 - x_17499*x_9413             external: the multiplier x_17499 = p
  a23618 = x_4432 - x_19964 - x_28730           external: L  = x_4432 - x_19964     TUNABLE
  a36659 = x_29854 - x_22665*x_1329             external: the multiplier x_22665 = p
  a36660 = 5113045*(x_7075*x_9118) - x_29854    external: K2 = 5113045*x_7075*x_9118 TUNABLE
  a36661 = x_31864 - x_28961*x_10903            external: the multiplier x_28961 = p
  a36662 = x_7075*x_8731                        external: J  = x_7075*x_8731        TUNABLE
  a36663 = x_31864                              external: NONE (purely private)
  a36664 = x_642 - x_28599*x_17325              external: the multiplier x_28599 = p

A configuration scan varies K1, L, K2, J.  It does NOT vary the four multipliers -- those are
all exactly p, which is a literal of the instance.  So the question is whether an admissible
boundary change exists with support inside {a23616, a23618, a36660, a36662}.
"""
import sys, json, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/tunable.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const0, cols = G.build_model(R, Pv, G.V0)
Eqs, rows0 = G.eq_system(R, Pv, const0, cols)
Bcol = {a: {e: G.EQCO[a].get(e, 0) for e in Eqs if G.EQCO[a].get(e, 0)} for a in sorted(const0)}
b0 = {e: rows0[e][1] for e in Eqs}

TUNABLE = [23616, 23618, 36660, 36662]
PMULT = [23617, 36659, 36661, 36664]
PRIVATE_ONLY = [36663]


def try_support(S):
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
    return (sol, msg)


say('TUNABLE by configuration : %s' % TUNABLE)
say('p-multiplier atoms       : %s  (external ingredient is the literal p)' % PMULT)
say('purely private           : %s' % PRIVATE_ONLY)

say('\n--- THE decisive test: boundary change confined to the configuration-tunable constants')
sol, msg = try_support(TUNABLE)
if sol is None:
    say('  UNSOLVABLE over Z: %s' % msg[:100])
    say('  => no configuration, of the 2,800 or of all 13,884, can clear this region.')
else:
    say('  SOLVABLE.  required change:')
    for a in TUNABLE:
        v = sol.get('d%d' % a, 0)
        say('     const(a%-6d) += %d bits %s' % (a, abs(v).bit_length(), '(zero)' if v == 0 else ''))

say('\n--- all subsets, smallest first: which supports admit an integral boundary change?')
works = []
allc = sorted(const0)
for k in range(1, len(allc) + 1):
    hit = []
    for S in itertools.combinations(allc, k):
        s, m = try_support(list(S))
        if s is not None:
            hit.append(S)
    say('  size %d: %d of %d supports work' % (k, len(hit), len(list(itertools.combinations(allc, k)))))
    for S in hit[:8]:
        tun = all(a in TUNABLE for a in S)
        say('      %s   %s' % (list(S), 'ALL TUNABLE' if tun else
                               'needs ' + str([a for a in S if a not in TUNABLE])))
    if hit:
        works = hit
        break

say('\n--- minimum number of p-multiplier atoms that must change')
if works:
    need = min(len([a for a in S if a in PMULT]) for S in works)
    say('  minimal supports of size %d; each needs at least %d p-multiplier atom(s)'
        % (len(works[0]), need))
say('DONE')
