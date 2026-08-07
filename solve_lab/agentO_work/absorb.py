"""Agent H's zero-collateral census lists NINE knobs
   {642,1329,8731,9118,9413,10903,17325,29854,31864}; my atom-level region model has eight
   {642,1329,9413,10903,17325,28730,29854,31864}.  The difference is x_8731 and x_9118, which
   my model holds fixed because they occur in atoms outside the region (1633/9081/40630 and
   1639/9083/36969/40630).  In H's frame those atoms are DEFINITIONS, so moving the knob just
   re-derives a downstream variable and the atom stays zero.

   To reproduce that inside my model I have to absorb those atoms into the region and check
   whether their downstream variables become private (so their residuals can be driven back
   to zero).  x_8731 and x_9118 enter the region's constants J and K2 directly, so if they
   become knobs the blocking congruences move.
"""
import sys, json, itertools, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/absorb.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


for u in (8731, 9118):
    say('x_%d occurs in atoms:' % u)
    for a in sorted(H.occ[u]):
        say('   a%-6d "%s"' % (a, H.atoms[a][:100]))
say('')
for u in (19892, 3349, 32010, 34310, 25297, 21279, 14865, 31861):
    say('x_%-6d occurs in %d atoms %s' % (u, len(H.occ[u]), sorted(H.occ[u])[:10]))
say('')

BASE = G.R0 + [23618]


def report(R, tag):
    R = sorted(set(R))
    Pv = G.private_vars(R)
    Eqs = sorted({e for x in R for e in G.EQCO[x]})
    say('%s: |R|=%d  private=%s  |E|=%d' % (tag, len(R), Pv, len(Eqs)))
    m = G.build_model(R, Pv, G.V0)
    if m is None:
        say('    NONLINEAR in a private var -> model not affine, skipped')
        return None
    const, cols = m
    Eqs, rows = G.eq_system(R, Pv, const, cols)
    t0 = time.time()
    sol, msg, _ = sparse.solve_sparse([rows[e][0] for e in Eqs], [rows[e][1] for e in Eqs],
                                      names=Eqs, verbose=False, maxcore=150,
                                      maxbits=10 ** 7, maxcorebits=10 ** 7)
    if sol is None:
        say('    FULL SYSTEM unsolvable over Z: %s  (%.0fs)' % (msg[:80], time.time() - t0))
        k, S, s2 = G.maxsat(Eqs, rows)
        say('    maxsat = %d of %d -> cost %d' % (k, len(Eqs), len(Eqs) - k))
        if len(Eqs) - k < 7:
            w = G.realise(R, Pv, s2)
            ff = E.eqfails(E.badatoms(w))
            say('    *** exact re-evaluation: %d failing -> score %d' % (len(ff), 39033 - len(ff)))
            if len(ff) < 7:
                json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0},
                          open('%s/absorb_%s_%d.json' % (OD, tag.replace('+', '_'), 39033 - len(ff)), 'w'))
                say('    *** WROTE improvement')
        return len(Eqs) - k
    w = G.realise(R, Pv, sol)
    ff = E.eqfails(E.badatoms(w))
    say('    INTEGRALLY SOLVABLE -> exact fails=%d score=%d  (%.0fs)'
        % (len(ff), 39033 - len(ff), time.time() - t0))
    if len(ff) < 7:
        json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0},
                  open('%s/absorb_%s_%d.json' % (OD, tag.replace('+', '_'), 39033 - len(ff)), 'w'))
        say('    *** WROTE improvement')
    return len(ff)


A8731 = sorted(set(H.occ[8731]) - set(BASE))
A9118 = sorted(set(H.occ[9118]) - set(BASE))
say('atoms to absorb for x_8731: %s' % A8731)
say('atoms to absorb for x_9118: %s' % A9118)
say('')
report(BASE, 'BASE')
report(BASE + A8731, 'BASE+free(x_8731)')
report(BASE + A9118, 'BASE+free(x_9118)')
report(BASE + A8731 + A9118, 'BASE+free(x_8731,x_9118)')
say('DONE')
