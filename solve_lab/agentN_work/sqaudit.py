"""Audit flagged by agent T: eq8680 is S^4, not S^2, so `inner` strips only one level.

atom 37887 parses as `S * S` with the two factors textually identical, so av[37887] = S^2 and
eq_terms[8680] = (1, sq=True, [(1, 37887)]) means eq8680 = (av[37887])^2 = S^4.  `optN.inner`
returns the sum over atom values, i.e. S^2 -- a QUADRATIC.  Zero locus is unaffected
(S^2 = 0 <=> S = 0), but any LINEAR model of that row is linearising a quadratic.

Two jobs here:
  1. find every atom in every row of my models that is a top-level square, so the blast radius is
     known rather than assumed;
  2. test affineness NUMERICALLY on every row of every model -- probe each knob at steps 1, 2, 3
     and require value(t) - value(0) to be exactly t * (value(1) - value(0)).  That catches any
     nonlinearity, not just the squares I know how to spot.
Then re-price the affected placements with the row stripped twice.
"""
import ast, json, itertools, re, sys, os
from collections import defaultdict
import ev, model
import optN
from optN import make, build, inner, WIT, POOL, atom_eqs
from widen import wide_knobs, build_wide
import zsolve

HERE = os.path.dirname(os.path.abspath(__file__))
VAR_RE = re.compile(r'x_(\d+)')
d = model.get()
atom_src = d['atom_src']
eq_terms = d['eq_terms']


def square_base(a):
    """if atom a is a top-level square, return the source of its base, else None"""
    t = ast.parse(atom_src[a], mode='eval').body
    if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Pow):
        try:
            if ast.literal_eval(t.right) == 2:
                return ast.unparse(t.left)
        except Exception:
            return None
    if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Mult):
        L, R = ast.unparse(t.left), ast.unparse(t.right)
        if L == R:
            return L
    return None


SQBASE = {}


def sqbase_code(a):
    if a not in SQBASE:
        b = square_base(a)
        SQBASE[a] = compile(VAR_RE.sub(r'v[\1]', b), '<s>', 'eval') if b else None
    return SQBASE[a]


def linear_core(st, e):
    """value of the correct LINEAR core of equation e, or None if the row is genuinely nonlinear.
    eq = m*(s^2 if sq else s), s = sum c*av[a];  zero <=> s = 0.
    If s is a single atom term and that atom is a square, strip again (recursively)."""
    m, sq, tl = eq_terms[e]
    terms = [(c, a) for c, a in tl]
    squares = [a for c, a in terms if square_base(a) is not None]
    if not squares:
        return inner(st, e), 0
    if len(terms) == 1:
        a = terms[0][1]
        depth = 1
        code = sqbase_code(a)
        val = eval(code, st.ns)
        # the base may itself be a square
        while True:
            src = square_base(a)
            b2 = None
            try:
                t2 = ast.parse(src, mode='eval').body
                if isinstance(t2, ast.BinOp) and isinstance(t2.op, ast.Mult) and \
                        ast.unparse(t2.left) == ast.unparse(t2.right):
                    b2 = ast.unparse(t2.left)
            except Exception:
                b2 = None
            if b2 is None:
                break
            code = compile(VAR_RE.sub(r'v[\1]', b2), '<s>', 'eval')
            val = eval(code, st.ns)
            depth += 1
            src = b2
        return val, depth
    return None, -1        # several atoms, at least one a square: no clean linear core


def affine_audit(D, tag):
    st = make(list(D))
    b0 = build(st)
    Rl = b0['R']
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows = list(Rl) + list(outside)
    print('\n=== %s ===  |R|=%d  outside=%d  wide knobs=%d' % (tag, len(Rl), len(outside), len(knobs)),
          flush=True)
    sqrows = []
    for e in rows:
        m, sq, tl = eq_terms[e]
        sqa = [a for c, a in tl if square_base(a) is not None]
        if sqa:
            sqrows.append((e, len(tl), sqa, e in Rl))
    print('  rows containing a top-level SQUARE atom: %d' % len(sqrows), flush=True)
    for e, nt, sqa, inR in sqrows:
        m, sq, tl = eq_terms[e]
        lc, dep = linear_core(st, e)
        print('     eq %-6d %s  sq=%-5s atom terms=%-3d square atoms=%s  -> %s'
              % (e, 'REGION' if inR else 'outside', sq, nt, sqa,
                 ('strip depth %d, core %s' % (dep, 'ZERO' if lc == 0 else 'nonzero'))
                 if lc is not None else 'NO CLEAN LINEAR CORE'), flush=True)
    # numeric affineness of every row against every knob
    bad = []
    base = [inner(st, e) for e in rows]
    for j, Y in enumerate(knobs):
        st1 = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        st2 = st.clone().set_free({Y: st.fv.get(Y, 0) + 2})
        st3 = st.clone().set_free({Y: st.fv.get(Y, 0) + 3})
        for i, e in enumerate(rows):
            d1 = inner(st1, e) - base[i]
            d2 = inner(st2, e) - base[i]
            d3 = inner(st3, e) - base[i]
            if d2 != 2 * d1 or d3 != 3 * d1:
                bad.append((e, Y, e in Rl))
    print('  NUMERIC affineness: %d (row, knob) pairs are NOT affine' % len(bad), flush=True)
    br = sorted(set(e for e, y, r in bad))
    print('  rows involved: %s' % br, flush=True)
    print('  of those, in the REGION: %s' % sorted(set(e for e, y, r in bad if r)), flush=True)
    return Rl, rows, bad, br


def reprice_stripped(D, tag):
    """re-price with row 8680 replaced by its true linear core S"""
    st = make(list(D))
    b0 = build(st)
    Rl, M, b, n = b0['R'], b0['M'], b0['b'], b0['n']
    knobs = b0['knobs']
    if 8680 not in Rl:
        print('  %s: 8680 not in region, nothing to strip' % tag, flush=True)
        return None
    i = Rl.index(8680)
    core0, dep = linear_core(st, 8680)
    newb = list(b)
    newb[i] = core0
    newM = [list(r) for r in M]
    for j, Y in enumerate(knobs):
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        c1, _ = linear_core(h, 8680)
        newM[i][j] = c1 - core0
    # affineness of the stripped row
    ok = True
    for j, Y in enumerate(knobs):
        h2 = st.clone().set_free({Y: st.fv.get(Y, 0) + 2})
        c2, _ = linear_core(h2, 8680)
        if c2 - core0 != 2 * newM[i][j]:
            ok = False
    opt_old, _, _, _ = zsolve.max_zero_rows(M, b, n, len(Rl))
    opt_new, rows_new, exh, _ = zsolve.max_zero_rows(newM, newb, n, len(Rl))
    print('  %-26s strip depth=%d  row-8680 affine after stripping: %-5s  '
          'OPT before=%d after=%d  failing=%d score=%d exh=%s'
          % (tag, dep, ok, opt_old, opt_new, len(Rl) - opt_new, 39033 - (len(Rl) - opt_new), exh),
          flush=True)
    return dict(D=list(D), strip_depth=dep, affine=ok, opt_old=opt_old, opt_new=opt_new,
                failing=len(Rl) - opt_new, score=39033 - (len(Rl) - opt_new), exh=str(exh))


if __name__ == '__main__':
    print('atom 37887: square_base is %s' % ('PRESENT (atom = S*S)' if square_base(37887) else 'absent'))
    print('eq_terms[8680] = %s' % (eq_terms[8680],))
    print('=> eq8680 = (av[37887])^2 = S^4 ; optN.inner returns av[37887] = S^2 (QUADRATIC)')

    affine_audit(WIT, 'WITNESS %s (|R|=12, 8680 absent)' % WIT)
    affine_audit([642], 'D=[642] (|R|=13, 8680 present)')

    print('\n=== re-price every |R|=13 state of the detach closure with 8680 stripped twice ===')
    KEY = [642, 28730, 29854, 31864]
    out = []
    for r in range(5):
        for D in itertools.combinations(KEY, r):
            if 28730 in D:
                continue          # |R|=12, 8680 not in the region
            res = reprice_stripped(list(D), 'D=%s' % (list(D),))
            if res:
                out.append(res)
    json.dump(out, open(os.path.join(HERE, 'runs', 'sqaudit.json'), 'w'), indent=1)
    if out:
        print('\nbest |R|=13 score after the correction: %d (was 39025)'
              % max(r['score'] for r in out))
