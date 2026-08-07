"""Single-atom region growth, restricted to atoms that actually FREE a new private variable
   (adding an atom that frees nothing can only add equations, never satisfy more)."""
import sys, json, time, signal
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/growB.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


R0 = G.R0
P0 = set(G.private_vars(R0))
E0 = sorted({e for x in R0 for e in G.EQCO[x]})
say('baseline: P=%d |E|=%d cost 7' % (len(P0), len(E0)))

cand = set()
for a in R0:
    for u in H.avars[a]:
        cand |= set(H.occ[u])
cand -= set(R0)
useful = []
for a in sorted(cand):
    P = set(G.private_vars(R0 + [a]))
    Eq = sorted({e for x in R0 + [a] for e in G.EQCO[x]})
    dP = P - P0
    dE = len(Eq) - len(E0)
    if dP:
        useful.append((a, sorted(dP), dE, len(Eq)))
say('%d adjacent atoms, %d of them free a new private var' % (len(cand), len(useful)))
for a, dP, dE, ne in useful:
    say('   a%d: frees %s, +%d equations (|E|=%d)  "%s"' % (a, dP, dE, ne, H.atoms[a][:60]))


class TO(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(TO()))
res = []
for a, dP, dE, ne in sorted(useful, key=lambda t: t[2]):
    if ne > 22:
        say('  +a%d: SKIP |E|=%d' % (a, ne))
        continue
    t0 = time.time()
    signal.alarm(600)
    try:
        rr = G.evaluate(R0 + [a], verbose=False)
    except TO:
        signal.alarm(0)
        say('  +a%d: TIMEOUT |E|=%d' % (a, ne))
        continue
    except Exception as ex:
        signal.alarm(0)
        say('  +a%d: ERR %s %s' % (a, type(ex).__name__, str(ex)[:60]))
        continue
    signal.alarm(0)
    if rr is None:
        say('  +a%d: NONLINEAR' % a)
        continue
    cost, ne2, k, S, sol, PP, RR = rr
    res.append((cost, a, ne2, k, len(PP)))
    say('  +a%d: P=%d |E|=%d maxsat=%d COST=%d (%.0fs)%s'
        % (a, len(PP), ne2, k, cost, time.time() - t0, ' ***BETTER***' if cost < 7 else ''))
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
json.dump(res, open(OD + '/grow2_results.json', 'w'))
say('DONE')
