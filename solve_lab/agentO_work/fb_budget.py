"""The one open door: a deliberately budgeted multi-atom compensation inside T.

T's 20 atoms all live in 10-18 equations (no free compensator), and the equations they disturb
are the region's own.  The 26 carriers of T = carriers of a37887 are therefore the entire
compensation channel, and they were all in the 34-knob set.  So what remains is purely the
BUDGET.  To beat 39,026 with j of the seven bought and b satisfied rows broken we need b < j.

  j=1, b=0 : already tested exhaustively -> none
  j=2, b<=1: 21 pairs x (no drop + each of the ~169 satisfied rows) -- COMPLETE, run here
  j=3, b<=2: not exhaustible; run a greedy upper bound and STATE it as a budget, not a proof
"""
import sys, os, json, re, itertools, time
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_budget.log', 'w', buffering=1)


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
say('knobs %d  rows %d  satisfied %d  failing %d' % (len(KNOB), len(names), len(SAT), len(FAIL)))
say('T = 0 is represented by the "S" row; breaking it = breaking eq8680.')


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
        out = OD + '/fbbudget_%d.json' % st.score()
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        say('   *** WROTE %s' % out)
        return st
    return None


say('\n=== j=2, b<=1 : COMPLETE test (buy 2 of the seven, break at most 1 arbitrary row)')
t0 = time.time()
pairs_feasible, nsolve, win = [], 0, None
for P in itertools.combinations(FAIL, 2):
    if solve(list(P)) is None:
        nsolve += 1
        continue                      # P alone infeasible => every superset infeasible
    nsolve += 1
    pairs_feasible.append(P)
say('pairs feasible on their own: %d of 21  -> %s' % (len(pairs_feasible), pairs_feasible))
for P in pairs_feasible:
    s = solve(SAT + list(P))
    nsolve += 1
    if s is not None:
        st = price(s, 'buy %s, break nothing' % list(P))
        if st:
            win = st
            break
    for r in SAT:
        keep = [e for e in SAT if e != r]
        s = solve(keep + list(P))
        nsolve += 1
        if s is not None:
            st = price(s, 'buy %s, break eq%s' % (list(P), r))
            if st:
                win = st
                break
    if win:
        break
say('j=2 complete test: %d solves, %.0fs, improvement found: %s'
    % (nsolve, time.time() - t0, win is not None))

if win is None:
    say('\n=> COMPLETE for j=2: no way to buy two of the seven while breaking at most one row.')
    say('   Combined with j=1,b=0 (exhaustive), the trade cannot be leveraged at budget 2.')

say('\n=== j>=3 : greedy upper bound (a BUDGET, not a proof)')
t0 = time.time()
for k in (3, 4, 5):
    for P in itertools.combinations(FAIL, k):
        if solve(list(P)) is None:
            continue
        keep = list(P)
        for r in SAT:
            cand = keep + [r]
            if solve(cand) is not None:
                keep = cand
        drops = len(SAT) - (len(keep) - k)
        say('   buy %s: greedy keeps %d of %d satisfied rows -> drops %d (need < %d)'
            % (list(P), len(keep) - k, len(SAT), drops, k))
        if drops < k:
            s = solve(keep)
            if s is not None:
                price(s, 'greedy buy %s' % list(P))
        break
    if time.time() - t0 > 1200:
        say('   (budget reached at k=%d)' % k)
        break
say('\nBEST: 39026 unless a WROTE line appears above')
say('DONE')
