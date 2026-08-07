"""j=3, b<=2 : exhaustive.  Buy three of the seven, break at most two rows -> 39,027 or better.

The greedy pass found that buying [12231,12270,12350] drops exactly 3 satisfied rows -- net
zero.  Greedy is only an UPPER bound on the number of drops, so the true minimum may be 2.
This enumerates it properly: for every triple, b=0 (1 solve), b=1 (168), b=2 (all C(168,2)).
Early exit on the first improvement.  Progress is logged per triple so the budget actually
tested is on the record even if the cap is hit.
"""
import sys, os, json, re, itertools, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_j3.log', 'w', buffering=1)
CAP = float(os.environ.get('OCAP', '2700'))


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
EQT = frameB.eq_terms
rows, rhs, names = {}, {}, []
for e in sorted({x for a in TOUCHED for x in fr.eq_of[a]}):
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
say('knobs %d  rows %d  satisfied %d  failing %d' % (len(KNOB), len(names), len(SAT), len(FAIL)))
say('target: buy 3, break <= 2  ->  39,027.   cap %.0fs' % CAP)


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                  names=[str(x) for x in S], verbose=False, maxcore=600,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
    return s


def price(s, tag):
    ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in KNOB if s.get(u, 0)}
    st = st0.clone().set_free(ch)
    say('   %s -> EXACT score %d  failing %s' % (tag, st.score(), sorted(st.fails)))
    if st.score() > 39026:
        out = OD + '/fbj3_%d.json' % st.score()
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('   *** WROTE %s' % out)
        return st
    return None


t0 = time.time()
nsolve = 0
win = None
done_triples = []
for P in itertools.combinations(FAIL, 3):
    if win or time.time() - t0 > CAP:
        break
    tp = time.time()
    if solve(list(P)) is None:
        nsolve += 1
        say(' triple %s: infeasible on its own -- skipped' % list(P))
        done_triples.append((P, 'infeasible-alone'))
        continue
    nsolve += 1
    # b = 0
    s = solve(SAT + list(P))
    nsolve += 1
    if s is not None:
        win = price(s, 'buy %s break 0' % list(P))
        done_triples.append((P, 'b=0 feasible'))
        continue
    # b = 1
    hit = False
    for r in SAT:
        s = solve([e for e in SAT if e != r] + list(P))
        nsolve += 1
        if s is not None:
            win = price(s, 'buy %s break eq%s' % (list(P), r))
            hit = True
            break
    if hit:
        done_triples.append((P, 'b=1 feasible'))
        continue
    # b = 2
    for r1, r2 in itertools.combinations(SAT, 2):
        if time.time() - t0 > CAP:
            break
        s = solve([e for e in SAT if e != r1 and e != r2] + list(P))
        nsolve += 1
        if s is not None:
            win = price(s, 'buy %s break eq%s,eq%s' % (list(P), r1, r2))
            hit = True
            break
    done_triples.append((P, 'b<=2 exhausted' if not hit else 'b=2 feasible'))
    say(' triple %s: b<=2 %s  (%d solves so far, %.0fs)'
        % (list(P), 'FEASIBLE' if hit else 'exhausted, none', nsolve, time.time() - tp))
    if hit:
        break

say('\ntriples completed at b<=2: %d of %d' % (len([1 for _, s in done_triples if 'exhaust' in s or 'infeas' in s]),
                                               len(list(itertools.combinations(FAIL, 3)))))
for P, st in done_triples:
    say('   %s : %s' % (list(P), st))
say('total solves %d, elapsed %.0fs, improvement: %s' % (nsolve, time.time() - t0, win is not None))
say('DONE')
