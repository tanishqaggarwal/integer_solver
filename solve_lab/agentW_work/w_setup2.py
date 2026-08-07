"""W round 2: SAME setup but with the EXTENDED knob set K+ = 38.
Round-1 K was 'free inputs reaching a NONZERO atom'; to repair a broken equation you may
instead move a DIFFERENT atom of that equation, and those movers were treated as constants.
K+ adds them.  (Original header follows.)

Rebuild O's frame-B budget setup verbatim (from fb_j3.py), as an importable module.

Produces: fr, st0, KNOB(34), rows/rhs over the TOUCHED equations, names, SAT(168+S), FAIL(7).
Only change vs O: paths point at agentW_work, and nothing is written into agentH_work.
"""
import sys, os, json, re, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.set_int_max_str_digits(20_000_000)
import frameB
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')
import sparse

VAR_RE = re.compile(r'x_(\d+)')
WITNESS = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'

fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open(WITNESS))
vw = [0] * frameB.NV
for k, val in W.items():
    vw[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv0 = {u: vw[u] for u in fr.free if vw[u] != 0}
st0 = frameB.State(fr, fv0)
assert st0.score() == 39026, st0.score()

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
KNOB = sorted(U | set(fr.SUPV.get(A, [])) | {u for e in (12231,12270,12350,14584,18673,22044,29125) for c,a in frameB.eq_terms[e][2] for u in fr.SUPV.get(a, [])})

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
rows, rhs, names, dropped = {}, {}, [], []
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
        dropped.append(e)
        continue
    rows[e] = {k: x for k, x in r.items() if x}
    rhs[e] = -const
    names.append(e)
rows['S'] = dict(Scol)
rhs['S'] = -S0
names.append('S')
SAT = [e for e in names if e == 'S' or e not in st0.fails]
FAIL = [e for e in names if e != 'S' and e in st0.fails]


def solve(S):
    s, m, _ = sparse.solve_sparse([rows[e] for e in S], [rhs[e] for e in S],
                                  names=[str(x) for x in S], verbose=False, maxcore=600,
                                  maxbits=10 ** 7, maxcorebits=10 ** 7)
    return s


def price(s, tag, log=None, tagfile='w'):
    ch = {u: fv0.get(u, 0) + s.get(u, 0) for u in KNOB if s.get(u, 0)}
    st = st0.clone().set_free(ch)
    msg = '   %s -> EXACT score %d  failing %s' % (tag, st.score(), sorted(st.fails))
    print(msg, file=log, flush=True) if log else print(msg, flush=True)
    if st.score() > 39026:
        out = os.path.join(HERE, 'w_%s_%d.json' % (tagfile, st.score()))
        json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
                  open(out, 'w'))
        print('   *** WROTE %s' % out, file=log, flush=True) if log else None
        return st
    return None


if __name__ == '__main__':
    print('knobs %d  rows %d  satisfied %d  failing %d  dropped-nonlinear %d'
          % (len(KNOB), len(names), len(SAT), len(FAIL), len(dropped)))
    print('KNOB =', KNOB)
    print('FAIL =', FAIL)
    print('dropped =', dropped)
    print('S row support =', len(Scol), 'S0 =', S0)
    print('reachable checks =', len(TOUCHED), 'reachable eqs =',
          len({x for a in TOUCHED for x in fr.eq_of[a]}))
