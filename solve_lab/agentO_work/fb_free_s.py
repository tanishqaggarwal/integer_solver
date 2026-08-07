"""Is there a knob that moves S independently of x_4432 and x_28730?

S = 0 is the binding constraint, and dS/dx_4432 = +1, dS/dx_28730 = -1 among the current 12
knobs.  If some other free input moves S at acceptable collateral, the L direction is restored
and delta0 becomes reachable.  Enumerate every free input supporting a37887, price it, and
re-run the exact maxsat with the cheap ones added.
"""
import sys, os, json, re, itertools, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_free_s.log', 'w', buffering=1)


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

BASEK = [642, 1329, 9413, 10903, 17325, 28730, 29854, 31864, 7068, 4432, 8731, 9118]
base = {a: st0.av[a] for a in fr.checks}
REGCHK = set()
for u in BASEK:
    o = fv0.get(u, 0)
    s1 = st0.clone().set_free({u: o + 1})
    REGCHK |= {a for a in fr.checks if s1.av[a] != base[a]}
say('checks reached by the 12 knobs: %s' % sorted(REGCHK))

sup = fr.SUPV.get(A, [])
say('free inputs supporting a%d: %d' % (A, len(sup)))
cand = []
for u in sup:
    o = fv0.get(u, 0)
    s1 = st0.clone().set_free({u: o + 1})
    ds = Sval(s1) - S0
    if ds == 0:
        continue
    touched = {a for a in fr.checks if s1.av[a] != base[a]}
    out = touched - REGCHK
    eqs = set()
    for a in out:
        eqs |= set(fr.eq_of[a])
    cand.append((len(eqs), len(out), u, ds))
cand.sort()
say('\nfree inputs that MOVE S, by outside-equation collateral:')
for n_eq, n_ck, u, ds in cand:
    say('   x_%-6d dS=%s  outside-checks %3d  outside-equations %4d%s'
        % (u, ds if abs(ds) < 10 ** 9 else '%db' % abs(ds).bit_length(), n_ck, n_eq,
           '   <<< in the 12' if u in BASEK else ''))
newk = [u for n_eq, n_ck, u, ds in cand if u not in BASEK]
say('\ncandidates outside the current 12: %s' % newk)
say('cheapest such: %s' % ([(u, n) for n, c, u, d in cand if u not in BASEK][:5]))

EQT = frameB.eq_terms


def attempt(knobs, tag):
    cols, nonlin, Scol = {}, set(), {}
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
        e1, e2 = Sval(s1) - S0, Sval(s2) - S0
        if e2 != 2 * e1:
            say('   %s: S nonlinear in x_%d -> skip' % (tag, u))
            return None
        if e1:
            Scol[u] = e1
    TOUCHED = set()
    for u in knobs:
        TOUCHED |= set(cols[u])
    BADc = {a for u, a in nonlin}
    if BADc - {A}:
        say('   %s: extra nonlinear checks %s' % (tag, sorted(BADc - {A})))
    EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
    rows, rhs, names = {}, {}, []
    for e in EQS:
        m, sq, tl = EQT[e]
        r, const = {}, 0
        skip = False
        for c, a in tl:
            if a == A:
                continue
            if a in BADc:
                skip = True
                break
            const += c * base.get(a, 0)
            for u in knobs:
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
    say('   %s: knobs=%d rows=%d (sat %d, fail %d)' % (tag, len(knobs), len(names), len(SAT), len(FAIL)))

    def solve(S):
        s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                      names=[str(x) for x in S], verbose=False, maxcore=400,
                                      maxbits=10 ** 7, maxcorebits=10 ** 7)
        return s
    best = None
    for k in range(len(FAIL), 0, -1):
        hit = False
        for S in itertools.combinations(FAIL, k):
            s = solve(SAT + list(S))
            if s is None:
                continue
            ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in knobs if s.get(u, 0)}
            st = st0.clone().set_free(ch)
            say('   %s: bought %s -> score %d  failing %s' % (tag, list(S), st.score(), sorted(st.fails)))
            hit = True
            if best is None or st.score() > best[0]:
                best = (st.score(), st)
            break
        if hit:
            break
    if best is None:
        say('   %s: nothing buyable' % tag)
    elif best[0] > 39026:
        out = OD + '/fbS_%d.json' % best[0]
        json.dump({('x_%d' % i): str(best[1].v[i]) for i in range(frameB.NV) if best[1].v[i] != 0},
                  open(out, 'w'))
        say('   *** WROTE %s' % out)
    return best


say('\n--- add the cheapest S-movers one at a time')
t0 = time.time()
for n_eq, n_ck, u, ds in cand:
    if u in BASEK:
        continue
    attempt(BASEK + [u], 'x_%d (collateral %d eqs)' % (u, n_eq))
    if time.time() - t0 > 900:
        say('   (budget reached)')
        break
say('DONE')
