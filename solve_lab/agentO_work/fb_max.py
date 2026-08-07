"""Exact maxsat on the frame-B region system.

12 knobs reach 12 checks and 29 equations; the witness satisfies 22 of the 29 and fails 7, and
ALL 7 failures are reachable.  Rather than greedy, start from the set the witness already
satisfies -- delta = 0 witnesses its solvability -- and ask which of the 7 failures can be
added.  Any single addition is 39,027.  Then try every subset of the failures.
"""
import sys, os, json, time, itertools
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_max.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


import frameB
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')
import sparse

DET = [642, 28730, 29854, 31864]
fr = frameB.Frame(DET)
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * frameB.NV
for k, val in W.items():
    vw[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv0 = {u: vw[u] for u in fr.free if vw[u] != 0}
st0 = frameB.State(fr, fv0)
assert st0.score() == 39026

PRIV = [642, 1329, 9413, 10903, 17325, 28730, 29854, 31864]
CARR = [7068, 4432, 8731, 9118]
KNOB = PRIV + CARR
base = {a: st0.av[a] for a in fr.checks}
cols, nonlin = {}, set()
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
TOUCHED = set()
for u in KNOB:
    TOUCHED |= set(cols[u])
BAD = {a for u, a in nonlin}
say('nonlinear checks: %s' % sorted(BAD))
for u, a in sorted(nonlin):
    say('   knob x_%d is nonlinear on check %d : "%s"' % (u, a, frameB.atom_src[a][:80]))

EQT = frameB.eq_terms
EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
rows, rhs, names, skipped = {}, {}, [], []
for e in EQS:
    m, sq, tl = EQT[e]
    r, const, ok = {}, 0, True
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
    rows[e] = {k: x for k, x in r.items() if x}
    rhs[e] = -const
    names.append(e)
SATIS = [e for e in names if e not in st0.fails]
FAIL = [e for e in names if e in st0.fails]
say('rows %d  (witness satisfies %d, fails %d)  skipped %s'
    % (len(names), len(SATIS), len(FAIL), skipped))
say('skipped equations currently failing: %s' % [e for e in skipped if e in st0.fails])


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S], names=list(S),
                                  verbose=False, maxcore=400, maxbits=10 ** 7,
                                  maxcorebits=10 ** 7)
    return s


def price(sol):
    ch = {u: fv0.get(u, 0) + sol.get(u, 0) for u in KNOB if sol.get(u, 0)}
    st = st0.clone().set_free(ch)
    return st.score(), sorted(st.fails), st


say('\nsanity: the set the witness already satisfies must be solvable')
s = solve(SATIS)
say('   solvable: %s' % (s is not None))

say('\n--- add one failing equation to the witness-satisfied set (each success = 39,027)')
wins = []
for f in FAIL:
    s = solve(SATIS + [f])
    if s is None:
        say('   +eq%-6d unsolvable' % f)
        continue
    sc, ff, st = price(s)
    say('   +eq%-6d SOLVABLE -> exact score %d  failing %s' % (f, sc, ff))
    wins.append((sc, f, s))
    if sc > 39026:
        out = OD + '/fbmax_%d_%d.json' % (f, sc)
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('   *** WROTE %s' % out)

say('\n--- every subset of the 7 failures, largest first')
best = None
for k in range(len(FAIL), 0, -1):
    got = False
    for S in itertools.combinations(FAIL, k):
        s = solve(SATIS + list(S))
        if s is None:
            continue
        sc, ff, st = price(s)
        say('   +%s SOLVABLE -> exact score %d  failing %s' % (list(S), sc, ff))
        got = True
        if best is None or sc > best[0]:
            best = (sc, list(S), st)
        break
    if got:
        break
    say('   no solvable subset of size %d' % k)
if best:
    sc, S, st = best
    say('\nBEST: %d  (added %s)' % (sc, S))
    if sc > 39026:
        out = OD + '/fbmax_best_%d.json' % sc
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('*** WROTE %s' % out)
else:
    say('\nNo failing equation can be added: the witness is exactly optimal over these 12 knobs.')
say('DONE')
