"""Can the WHOLE region system be solved over Z?

For region R the only equations that can fail are E(R); everything else is untouched because
the private variables of R occur in no atom outside R.  So if the full system E(R) is
integrally solvable, the resulting assignment satisfies ALL 39,033 equations.

For R0 + {a23618} the system has rank 8 = #unknowns, the 5 dependent rows are exactly
consistent, so over Q the solution is unique — and exactly four divisibilities block it
(three by p, one by a 279-bit pivot).  Growing the region adds knobs; this script tests every
single- and double-atom growth for full integral solvability, and reports how many of the
blocking congruences survive.
"""
import sys, json, time, itertools, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fullreg.log', 'w', buffering=1)
P256 = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


def test(R, tag, verbose=True):
    P = G.private_vars(R)
    m = G.build_model(R, P, G.V0)
    if m is None:
        say('%s: NONLINEAR' % tag)
        return None
    const, cols = m
    Eqs, rows = G.eq_system(R, P, const, cols)
    t0 = time.time()
    sol, msg, _ = sparse.solve_sparse([rows[e][0] for e in Eqs], [rows[e][1] for e in Eqs],
                                      names=Eqs, verbose=False, maxcore=120,
                                      maxbits=10 ** 7, maxcorebits=10 ** 7)
    if sol is None:
        if verbose:
            say('%s: P=%d |E|=%d  FULL SYSTEM UNSOLVABLE over Z: %s  (%.0fs)'
                % (tag, len(P), len(Eqs), msg[:70], time.time() - t0))
        return ('unsat', len(P), len(Eqs), msg)
    w = G.realise(R, P, sol)
    ff = E.eqfails(E.badatoms(w))
    say('%s: P=%d |E|=%d  INTEGRALLY SOLVABLE -> exact fails=%d score=%d  (%.0fs)'
        % (tag, len(P), len(Eqs), len(ff), 39033 - len(ff), time.time() - t0))
    if len(ff) < 7:
        out = '%s/full_%s_%d.json' % (OD, tag.replace(' ', '').replace(',', '_'), 39033 - len(ff))
        json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0}, open(out, 'w'))
        say('    *** WROTE %s' % out)
    return ('sat', len(P), len(Eqs), len(ff))


BASE = G.R0 + [23618]
say('p is the 256-bit constant: %s...' % str(P256)[:24])
test(BASE, 'R0+a23618')

cands = [c['atom'] for c in json.load(open(OD + '/growcand.json')) if c['atom'] != 23618]
say('\n--- single growths on top of R0+a23618 (%d candidates)' % len(cands))
sat = []
for a in cands:
    r = test(BASE + [a], 'R0+a23618+a%d' % a)
    if r and r[0] == 'sat':
        sat.append((a, r[3]))
say('\nsolvable single growths: %s' % (sat,))

say('\n--- double growths')
best = None
t0 = time.time()
for a, b in itertools.combinations(cands, 2):
    r = test(BASE + [a, b], 'R0+a23618+a%d+a%d' % (a, b), verbose=False)
    if r and r[0] == 'sat':
        say('  PAIR (%d,%d) solvable, fails=%d' % (a, b, r[3]))
        if best is None or r[3] < best[2]:
            best = (a, b, r[3])
    if time.time() - t0 > 5400:
        say('  (double-growth budget exhausted)')
        break
say('best double: %s' % (best,))
say('DONE')
