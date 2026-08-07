"""Single-atom region growth around the 39,026 witness residual, with guards.
   Writes its own log (shell redirection into task files has proved unreliable here)."""
import sys, json, time, signal
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/growA.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


R0 = G.R0
r = G.evaluate(R0, verbose=False)
say('baseline: |E|=%d maxsat=%d COST=%d  S=%s' % (r[1], r[2], r[0], r[3]))
w = G.realise(r[6], r[5], r[4])
say('  exact:', len(E.eqfails(E.badatoms(w))), 'failing equations')

cand = set()
for a in R0:
    for u in H.avars[a]:
        cand |= set(H.occ[u])
cand -= set(R0)
say('%d adjacent atoms' % len(cand))


class TO(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(TO()))
res = []
for a in sorted(cand):
    R = R0 + [a]
    P = G.private_vars(R)
    Eqs = sorted({e for x in R for e in G.EQCO[x]})
    if len(P) > 60 or len(Eqs) > 30:
        say('  +a%d: SKIP P=%d |E|=%d' % (a, len(P), len(Eqs)))
        continue
    t0 = time.time()
    signal.alarm(120)
    try:
        rr = G.evaluate(R, verbose=False)
    except TO:
        signal.alarm(0)
        say('  +a%d: TIMEOUT P=%d |E|=%d' % (a, len(P), len(Eqs)))
        continue
    except Exception as ex:
        signal.alarm(0)
        say('  +a%d: ERR %s %s' % (a, type(ex).__name__, str(ex)[:60]))
        continue
    signal.alarm(0)
    if rr is None:
        say('  +a%d: NONLINEAR P=%d |E|=%d' % (a, len(P), len(Eqs)))
        continue
    cost, ne, k, S, sol, PP, RR = rr
    res.append((cost, a, ne, k, len(PP)))
    say('  +a%d: P=%d |E|=%d maxsat=%d COST=%d (%.0fs)%s'
        % (a, len(PP), ne, k, cost, time.time() - t0, ' ***BETTER***' if cost < 7 else ''))
    if cost < 7:
        w = G.realise(RR, PP, sol)
        ff = E.eqfails(E.badatoms(w))
        say('      EXACT %d failing -> score %d' % (len(ff), 39033 - len(ff)))
        if len(ff) < 7:
            json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0},
                      open('%s/grow_%d_%d.json' % (OD, a, 39033 - len(ff)), 'w'))
            say('      *** WROTE improvement')
res.sort()
say('BEST: %s' % (res[:15],))
json.dump(res, open(OD + '/grow1_results.json', 'w'))
say('DONE')
