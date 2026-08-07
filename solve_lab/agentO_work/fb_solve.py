"""Realize delta0 natively in frame B and measure the collateral EXACTLY.

My atom-level model could not express re-derivation, so it had to treat the four boundary
constants as free-floating.  Frame B re-derives, and it turns out both open carriers are
themselves free inputs there (x_7068 supports the K1 atom, x_4432 the L atom).  So the honest
computation is to redo the region solve natively in frame B over

    the 8 region-private variables  +  the 4 carriers  x_7068, x_4432, x_8731, x_9118

and to take as ROWS every equation that any of those twelve can actually reach -- not just my
13.  Then the solve either closes everything (score 39,033) or the shortfall IS the collateral,
measured rather than bounded.

agentH_work is imported READ-ONLY.
"""
import sys, os, json, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
sys.path.append(OD)
LOG = open(OD + '/runs/fb_solve.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


import frameB
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')  # sparse.py lives there
import sparse

DET = [642, 28730, 29854, 31864]
fr = frameB.Frame(DET)
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * frameB.NV
for k, val in W.items():
    vw[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv0 = {u: vw[u] for u in fr.free if vw[u] != 0}
st0 = frameB.State(fr, fv0)
say('witness in frame B: score %d  failing %s' % (st0.score(), sorted(st0.fails)))
assert st0.score() == 39026

PRIV = [642, 1329, 9413, 10903, 17325, 28730, 29854, 31864]
CARR = [7068, 4432, 8731, 9118]
KNOB = PRIV + CARR
for u in KNOB:
    assert u in set(fr.free), 'x_%d is not free in frame B' % u
say('knobs (all free in frame B): %s' % KNOB)

# --- exact columns: d(check atom)/d(knob), with a linearity check
base = {a: st0.av[a] for a in fr.checks}
cols, nonlin = {}, set()
t0 = time.time()
for u in KNOB:
    o = fv0.get(u, 0)
    s1 = st0.clone().set_free({u: o + 1})
    s2 = st0.clone().set_free({u: o + 2})
    c = {}
    for a in fr.checks:
        d1 = s1.av[a] - base[a]
        if d1:
            c[a] = d1
        if s2.av[a] - base[a] != 2 * d1:
            nonlin.add((u, a))
    cols[u] = c
say('columns probed in %.0fs; nonlinear (knob,check) pairs: %d' % (time.time() - t0, len(nonlin)))
for u in KNOB:
    nl = [a for a in cols[u] if (u, a) in nonlin]
    say('   x_%-6d touches %3d checks, %d of them nonlinearly' % (u, len(cols[u]), len(nl)))

TOUCHED = set()
for u in KNOB:
    TOUCHED |= set(cols[u])
BAD = {a for u, a in nonlin}
say('checks reachable by the knobs: %d   (nonlinear in some knob: %d)' % (len(TOUCHED), len(BAD)))
USE = sorted(TOUCHED - BAD)
say('checks usable as linear rows: %d' % len(USE))

# --- every equation any reachable check appears in
EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
say('equations reachable by the knobs: %d  (of which currently failing: %d)'
    % (len(EQS), len([e for e in EQS if e in st0.fails])))
say('witness failures inside this set: %s' % sorted(set(st0.fails) & set(EQS)))
say('witness failures OUTSIDE it (unreachable, permanent): %s' % sorted(set(st0.fails) - set(EQS)))

# --- rows: for equation e require sum_a c_{e,a} * av_a = 0
EQT = frameB.eq_terms
rows, rhs, names, skipped = [], [], [], []
for e in EQS:
    m, sq, tl = EQT[e]
    r = {}
    const = 0
    ok = True
    for c, a in tl:
        if a in BAD:
            ok = False
            break
        const += c * base.get(a, 0)
        for u in KNOB:
            k = cols[u].get(a)
            if k:
                r[u] = r.get(u, 0) + c * k
    if not ok:
        skipped.append(e)
        continue
    r = {k: x for k, x in r.items() if x}
    rows.append(r)
    rhs.append(-const)
    names.append(e)
say('rows built: %d   skipped (nonlinear check in the equation): %d' % (len(rows), len(skipped)))

t0 = time.time()
sol, msg, _ = sparse.solve_sparse(rows, rhs, names=names, verbose=False, maxcore=400,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
say('\nfull solve over Z: %s  (%.0fs)' % (msg[:90], time.time() - t0))
if sol is not None:
    ch = {u: fv0.get(u, 0) + sol.get(u, 0) for u in KNOB if sol.get(u, 0)}
    st = st0.clone().set_free(ch)
    say('APPLIED: score %d  failing %s' % (st.score(), sorted(st.fails)))
    if st.score() > 39026:
        out = '/home/user/integer_solver/solve_lab/agentO_work/fbsolve_%d.json' % st.score()
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('*** WROTE %s' % out)
else:
    say('not integrally solvable on the full row set; measuring the price by dropping rows')
    import itertools
    # greedy: keep as many equations as possible, prioritising the currently-failing ones
    order = [e for e in names if e in st0.fails] + [e for e in names if e not in st0.fails]
    idx = {e: i for i, e in enumerate(names)}
    keep, ksol = [], None
    for e in order:
        cand = keep + [e]
        s, m2, _ = sparse.solve_sparse([rows[idx[x]] for x in cand], [rhs[idx[x]] for x in cand],
                                       names=cand, verbose=False, maxcore=400,
                                       maxbits=10 ** 7, maxcorebits=10 ** 7)
        if s is not None:
            keep, ksol = cand, s
    say('greedy kept %d of %d equations' % (len(keep), len(names)))
    if ksol is not None:
        ch = {u: fv0.get(u, 0) + ksol.get(u, 0) for u in KNOB if ksol.get(u, 0)}
        st = st0.clone().set_free(ch)
        say('APPLIED greedy: score %d  failing %s' % (st.score(), sorted(st.fails)))
        if st.score() > 39026:
            out = '/home/user/integer_solver/solve_lab/agentO_work/fbsolve_%d.json' % st.score()
            json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                      open(out, 'w'))
            say('*** WROTE %s' % out)
say('DONE')
