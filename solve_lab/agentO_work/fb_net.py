"""Close the remaining gap: NET gain, not zero-collateral gain.

fb_U required every currently-satisfied row to stay satisfied.  A score above 39,026 does not
need that -- it needs (bought) > (broken).  The known tax is exactly one equation, eq8680
(the S row).  So the decisive test is: drop the S row -- i.e. allow eq8680 to break, paying 1 --
and try to buy TWO or more of the seven.  Any success is 39,027 or better.

Also tested: paying one arbitrary satisfied row instead, and paying two.
"""
import sys, os, json, re, itertools, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_net.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


import frameB
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')
import sparse

VAR_RE = re.compile(r'x_(\d+)')
fr = frameB.Frame([642, 28730, 29854, 31864])
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
SCODE = compile(VAR_RE.sub(r'v[\1]', src[:half]), '<S>', 'eval')
Sval = lambda st: eval(SCODE, {'v': st.v, '__builtins__': {}})
S0 = Sval(st0)

NZ = sorted(st0.nz())
U = set()
for a in NZ:
    U |= set(fr.SUPV.get(a, []))
KNOB = sorted(U | set(fr.SUPV.get(A, [])))
say('knobs: %d' % len(KNOB))

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
    if e2 == 2 * e1 and e1:
        Scol[u] = e1
TOUCHED = set()
for u in KNOB:
    TOUCHED |= set(cols[u])
BADc = {a for u, a in nonlin}
EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
EQT = frameB.eq_terms
rows, rhs, names = {}, {}, []
for e in EQS:
    m, sq, tl = EQT[e]
    r, const, skip = {}, 0, False
    for c, a in tl:
        if a == A:
            continue
        if a in BADc:
            skip = True
            break
        const += c * base.get(a, 0)
        for u in KNOB:
            k = cols[u].get(a)
            if k:
                r[u] = r.get(u, 0) + c * k
    if skip:
        continue
    rows[e] = {k: x for k, x in r.items() if x}
    rhs[e] = -const
    names.append(e)
rows['S'] = dict(Scol)
rhs['S'] = -S0
names.append('S')
SAT = [e for e in names if e == 'S' or e not in st0.fails]
FAIL = [e for e in names if e != 'S' and e in st0.fails]
say('rows %d  satisfied %d  failing %d' % (len(names), len(SAT), len(FAIL)))


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                  names=[str(x) for x in S], verbose=False, maxcore=600,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
    return s


def price(s):
    ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in KNOB if s.get(u, 0)}
    st = st0.clone().set_free(ch)
    return st.score(), sorted(st.fails), st


BEST = [39026, None]


def record(s, tag):
    sc, ff, st = price(s)
    say('   %s -> EXACT score %d  failing %s' % (tag, sc, ff))
    if sc > BEST[0]:
        BEST[0], BEST[1] = sc, st
        out = OD + '/fbnet_%d.json' % sc
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('   *** WROTE %s' % out)


say('\n=== TEST A: pay eq8680 (drop the S row), buy k >= 2 of the seven')
SAT_noS = [e for e in SAT if e != 'S']
for k in range(len(FAIL), 1, -1):
    hit = False
    for Sx in itertools.combinations(FAIL, k):
        s = solve(SAT_noS + list(Sx))
        if s is None:
            continue
        say('  bought %d: %s' % (k, list(Sx)))
        record(s, 'pay eq8680, buy %s' % list(Sx))
        hit = True
        break
    if hit:
        break
    say('  k=%d: none' % k)

say('\n=== TEST B: pay eq8680 + one other satisfied row, buy k >= 3')
t0 = time.time()
done = False
for k in (3, 4):
    if done:
        break
    for drop in SAT_noS:
        if time.time() - t0 > 1500:
            say('  (budget reached)')
            done = True
            break
        keep = [e for e in SAT_noS if e != drop]
        for Sx in itertools.combinations(FAIL, k):
            s = solve(keep + list(Sx))
            if s is None:
                continue
            say('  pay eq8680 + eq%s, bought %s' % (drop, list(Sx)))
            record(s, 'pay eq8680+eq%s, buy %s' % (drop, list(Sx)))
            done = True
            break
        if done:
            break
    if not done:
        say('  k=%d: none over any single extra payment' % k)

say('\nBEST OVERALL: %d' % BEST[0])
say('DONE')
