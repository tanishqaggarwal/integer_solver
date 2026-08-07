"""The whole barrier is now one atom: a37887, nonlinear in x_4432 and x_28730.

Every one of the seven failing equations can be bought, and every purchase costs eq8680 --
the single equation containing a37887, which had to be dropped from the linear row set.
Characterise that dependence exactly (degree, coefficients), then try to satisfy eq8680 too:
  (a) hold each nonlinear knob fixed and see whether the atom linearises;
  (b) if the dependence is quadratic, solve it exactly and feed the roots back.
"""
import sys, os, json, itertools
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_nl.log', 'w', buffering=1)


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
A = 37887
say('atom a%d full source:\n   %s\n' % (A, frameB.atom_src[A]))
say('a%d appears in equations: %s' % (A, fr.eq_of[A]))
say('a%d current value: %d bits' % (A, abs(st0.av[A]).bit_length()))

say('\n--- degree of a%d in each nonlinear knob' % A)
for u in (4432, 28730):
    o = fv0.get(u, 0)
    vals = []
    for t in range(5):
        s = st0.clone().set_free({u: o + t})
        vals.append(s.av[A])
    d1 = [vals[i + 1] - vals[i] for i in range(4)]
    d2 = [d1[i + 1] - d1[i] for i in range(3)]
    d3 = [d2[i + 1] - d2[i] for i in range(2)]
    deg = 0 if all(x == 0 for x in d1) else (1 if all(x == 0 for x in d2)
                                             else (2 if all(x == 0 for x in d3) else '>2'))
    say('   x_%-6d degree %s   2nd difference constant = %s'
        % (u, deg, d2[0] if deg == 2 else '-'))
    if deg == 2:
        a2 = d2[0] // 2
        a1 = d1[0] - a2
        a0 = vals[0]
        say('        a%d(x_%d = o+t) = %s*t^2 + %s*t + const   (leading coeff %d bits)'
            % (A, u, str(a2)[:20] + '..', str(a1)[:20] + '..', abs(a2).bit_length()))

PRIV = [642, 1329, 9413, 10903, 17325, 28730, 29854, 31864]
CARR = [7068, 4432, 8731, 9118]
EQT = frameB.eq_terms


def build(knobs):
    base = {a: st0.av[a] for a in fr.checks}
    cols, nonlin = {}, set()
    for u in knobs:
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
    for u in knobs:
        TOUCHED |= set(cols[u])
    BADc = {a for u, a in nonlin}
    EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
    rows, rhs, names, skipped = {}, {}, [], []
    for e in EQS:
        m, sq, tl = EQT[e]
        r, const, ok = {}, 0, True
        for c, a in tl:
            if a in BADc:
                ok = False
                break
            const += c * base.get(a, 0)
            for u in knobs:
                k = cols[u].get(a)
                if k:
                    r[u] = r.get(u, 0) + c * k
        if not ok:
            skipped.append(e)
            continue
        rows[e] = {k: x for k, x in r.items() if x}
        rhs[e] = -const
        names.append(e)
    return rows, rhs, names, skipped, BADc


def run(knobs, tag):
    rows, rhs, names, skipped, BADc = build(knobs)
    SAT = [e for e in names if e not in st0.fails]
    FAIL = [e for e in names if e in st0.fails]
    say('\n=== %s : knobs=%d rows=%d (sat %d, fail %d) skipped=%s nonlinear-checks=%s'
        % (tag, len(knobs), len(names), len(SAT), len(FAIL), skipped, sorted(BADc)))
    if 8680 not in names and 8680 in [e for a in BADc for e in fr.eq_of[a]]:
        say('    eq8680 still not linear -> cannot be constrained here')

    def solve(S):
        s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S], names=list(S),
                                      verbose=False, maxcore=400, maxbits=10 ** 7,
                                      maxcorebits=10 ** 7)
        return s
    best = None
    for f in FAIL:
        s = solve(SAT + [f])
        if s is None:
            continue
        ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in knobs if s.get(u, 0)}
        st = st0.clone().set_free(ch)
        say('    +eq%-6d -> score %d  failing %s' % (f, st.score(), sorted(st.fails)))
        if best is None or st.score() > best[0]:
            best = (st.score(), f, st)
    if best:
        say('    best here: %d' % best[0])
        if best[0] > 39026:
            out = OD + '/fbnl_%d.json' % best[0]
            json.dump({('x_%d' % i): str(best[2].v[i]) for i in range(frameB.NV) if best[2].v[i] != 0},
                      open(out, 'w'))
            say('    *** WROTE %s' % out)
    else:
        say('    nothing buyable')
    return best


run(PRIV + CARR, 'all 12 knobs (baseline)')
run([u for u in PRIV + CARR if u != 4432], 'x_4432 held fixed')
run([u for u in PRIV + CARR if u != 28730], 'x_28730 held fixed')
run([u for u in PRIV + CARR if u not in (4432, 28730)], 'both nonlinear knobs held fixed')
say('DONE')
