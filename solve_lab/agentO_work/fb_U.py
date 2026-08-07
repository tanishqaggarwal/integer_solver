"""Upgrade the seven-way trade from a knob-set measurement to a scoped THEOREM.

Facts now established:
  (F1) eq8680 has exactly one term, coefficient 1, squared-outer, on a37887 = S^2.
       So eq8680 = S^4, and S = 0 in EVERY satisfying assignment.  Unconditional.
  (F2) The 7 failing equations are functions of the witness's 8 nonzero check atoms alone.

(F2) means those 7 equations can only be affected by free inputs that reach one of those 8
atoms.  Call that set U.  Any assignment that agrees with the witness OUTSIDE U therefore
scores at most whatever the best assignment inside U scores.  So run the exact maxsat over the
whole of U (plus every carrier of S, so the S = 0 constraint is fully represented) and the
answer is optimal over that entire class -- not over a hand-picked 12.
"""
import sys, os, json, re, itertools, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_U.log', 'w', buffering=1)


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
say('witness nonzero check atoms (the residual): %s' % NZ)
FAILEQ = sorted(st0.fails)
say('witness failing equations: %s' % FAILEQ)

# sanity for (F2): the failing equations must involve only atoms in NZ (others are zero)
for e in FAILEQ:
    m, sq, tl = frameB.eq_terms[e]
    nzin = [a for c, a in tl if a in set(NZ)]
    say('   eq%-6d has %d terms, %d of them on nonzero atoms %s' % (e, len(tl), len(nzin), nzin))

U = set()
for a in NZ:
    U |= set(fr.SUPV.get(a, []))
say('\n|U| = free inputs reaching a nonzero region atom: %d' % len(U))
US = set(fr.SUPV.get(A, []))
say('carriers of S: %d' % len(US))
KNOB = sorted(U | US)
say('knob set = U + carriers(S): %d free inputs' % len(KNOB))

base = {a: st0.av[a] for a in fr.checks}
cols, nonlin, Scol = {}, set(), {}
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
    e1, e2 = Sval(s1) - S0, Sval(s2) - S0
    if e2 == 2 * e1 and e1:
        Scol[u] = e1
    elif e2 != 2 * e1:
        say('   S nonlinear in x_%d -- excluded from the S row' % u)
say('columns probed in %.0fs' % (time.time() - t0))

TOUCHED = set()
for u in KNOB:
    TOUCHED |= set(cols[u])
BADc = {a for u, a in nonlin}
say('checks reachable: %d   nonlinear: %d %s' % (len(TOUCHED), len(BADc), sorted(BADc)[:10]))
EQS = sorted({e for a in TOUCHED for e in fr.eq_of[a]})
say('equations reachable: %d   (witness failures inside: %d of %d)'
    % (len(EQS), len(set(FAILEQ) & set(EQS)), len(FAILEQ)))
miss = set(FAILEQ) - set(EQS)
say('failures NOT reachable (permanently failing for this class): %s' % sorted(miss))

EQT = frameB.eq_terms
rows, rhs, names, dropped = {}, {}, [], []
for e in EQS:
    m, sq, tl = EQT[e]
    r, const, skip = {}, 0, False
    for c, a in tl:
        if a == A:
            continue                    # forced to 0 by the S row
        if a in BADc:
            skip = True
            break
        const += c * base.get(a, 0)
        for u in KNOB:
            k = cols[u].get(a)
            if k:
                r[u] = r.get(u, 0) + c * k
    if skip:
        dropped.append(e)
        continue
    rows[e] = {k: x for k, x in r.items() if x}
    rhs[e] = -const
    names.append(e)
rows['S'] = dict(Scol)
rhs['S'] = -S0
names.append('S')
say('rows: %d   dropped (nonlinear atom): %d %s' % (len(names), len(dropped), dropped[:10]))
say('dropped rows that currently FAIL: %s' % [e for e in dropped if e in st0.fails])

SAT = [e for e in names if e == 'S' or e not in st0.fails]
FAIL = [e for e in names if e != 'S' and e in st0.fails]
say('witness satisfies %d rows, fails %d: %s' % (len(SAT), len(FAIL), FAIL))


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                  names=[str(x) for x in S], verbose=False, maxcore=600,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
    return s


say('\nsanity: witness-satisfied set solvable: %s' % (solve(SAT) is not None))
say('\n--- can ANY failing equation be bought over the whole of U + carriers(S)?')
best = None
for k in range(len(FAIL), 0, -1):
    hit = False
    for Sx in itertools.combinations(FAIL, k):
        s = solve(SAT + list(Sx))
        if s is None:
            continue
        ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in KNOB if s.get(u, 0)}
        st = st0.clone().set_free(ch)
        say('   bought %s -> score %d  failing %s' % (list(Sx), st.score(), sorted(st.fails)))
        hit = True
        if best is None or st.score() > best[0]:
            best = (st.score(), st)
        break
    if hit:
        break
    say('   no solvable subset of size %d' % k)

if best and best[0] > 39026:
    out = OD + '/fbU_%d.json' % best[0]
    json.dump({('x_%d' % i): str(best[1].v[i]) for i in range(frameB.NV) if best[1].v[i] != 0},
              open(out, 'w'))
    say('*** WROTE %s' % out)
elif best:
    say('\nbest achievable: %d (no improvement)' % best[0])
else:
    say('\nNOTHING buyable over the entire class.')
say('DONE')
