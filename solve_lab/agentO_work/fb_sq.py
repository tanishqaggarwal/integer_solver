"""a37887 is a PERFECT SQUARE: its source is literally (S) * (S), its value at the witness is 0,
and the second difference in each nonlinear knob is exactly 2, so S is linear with coefficient
+/-1 in each.  Therefore "a37887 = 0" is not a quadratic constraint -- it is the LINEAR
constraint S = 0, which can be added as an ordinary row.  That is why eq8680 had to be dropped,
and why every repair cost it.

Add S = 0 as a row and redo the exact maxsat.
"""
import sys, os, json, re, itertools
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_sq.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


import frameB
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')
import sparse

VAR_RE = re.compile(r'x_(\d+)')
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
src = frameB.atom_src[A]
half = (len(src) - 3) // 2
assert src[half:half + 3] == ' * ' and src[:half] == src[half + 3:], 'not a clean square'
Ssrc = src[:half]
say('a%d = (S)*(S) confirmed; S has %d chars' % (A, len(Ssrc)))
SCODE = compile(VAR_RE.sub(r'v[\1]', Ssrc), '<S>', 'eval')


def Sval(st):
    return eval(SCODE, {'v': st.v, '__builtins__': {}})


S0 = Sval(st0)
say('S at the witness = %s   (a%d = %s)' % (S0, A, st0.av[A]))
assert S0 * S0 == st0.av[A]

PRIV = [642, 1329, 9413, 10903, 17325, 28730, 29854, 31864]
CARR = [7068, 4432, 8731, 9118]
KNOB = PRIV + CARR
EQT = frameB.eq_terms
base = {a: st0.av[a] for a in fr.checks}

cols, nonlin, Scol = {}, set(), {}
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
    e1, e2 = Sval(s1) - S0, Sval(s2) - S0
    assert e2 == 2 * e1, 'S is not linear in x_%d' % u
    if e1:
        Scol[u] = e1
say('S is linear in every knob; dS/dknob = %s' % {('x_%d' % u): c for u, c in Scol.items()})

TOUCHED = set()
for u in KNOB:
    TOUCHED |= set(cols[u])
BADc = {a for u, a in nonlin}
say('nonlinear checks: %s   (handled analytically via S)' % sorted(BADc))
EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})

rows, rhs, names = {}, {}, []
for e in EQS:
    m, sq, tl = EQT[e]
    r, const, uses_sq = {}, 0, False
    for c, a in tl:
        if a in BADc:
            uses_sq = True
            continue          # a37887 is forced to 0 by the S row below
        const += c * base.get(a, 0)
        for u in KNOB:
            k = cols[u].get(a)
            if k:
                r[u] = r.get(u, 0) + c * k
    rows[e] = {k: x for k, x in r.items() if x}
    rhs[e] = -const
    names.append(e)
    if uses_sq:
        say('  eq%d contains a%d; with S=0 its remaining row is %s (rhs %s)'
            % (e, A, rows[e] if rows[e] else 'EMPTY', rhs[e]))
# the S row
rows['S'] = dict(Scol)
rhs['S'] = -S0
names.append('S')
say('rows: %d (including the S row)' % len(names))

SAT = [e for e in names if e == 'S' or e not in st0.fails]
FAIL = [e for e in names if e != 'S' and e in st0.fails]
say('witness satisfies %d rows, fails %d' % (len(SAT), len(FAIL)))


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                  names=[str(x) for x in S], verbose=False, maxcore=400,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
    return s


def price(sol):
    ch = {u: fv0.get(u, 0) + sol.get(u, 0) for u in KNOB if sol.get(u, 0)}
    st = st0.clone().set_free(ch)
    return st.score(), sorted(st.fails), st


say('\nsanity: witness-satisfied set (with S row) solvable: %s' % (solve(SAT) is not None))
say('\n--- buy failing equations, keeping every satisfied row AND S = 0')
best = None
for k in range(len(FAIL), 0, -1):
    found = False
    for S in itertools.combinations(FAIL, k):
        s = solve(SAT + list(S))
        if s is None:
            continue
        sc, ff, st = price(s)
        say('   + %s  -> score %d  failing %s' % (list(S), sc, ff))
        found = True
        if best is None or sc > best[0]:
            best = (sc, list(S), st)
        break
    if found:
        break
    say('   no solvable subset of size %d' % k)

if best:
    sc, S, st = best
    say('\nBEST %d (bought %s)' % (sc, S))
    if sc > 39026:
        out = OD + '/fbsq_%d.json' % sc
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('*** WROTE %s' % out)
else:
    say('\nnothing buyable even with S handled analytically')
say('DONE')
