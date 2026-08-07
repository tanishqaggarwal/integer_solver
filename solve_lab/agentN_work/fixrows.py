"""Redo every analysis whose rows were mis-modelled, with the square rows stripped twice.

The numeric affineness audit (sqaudit.py) found EXACTLY two non-affine rows anywhere in my models:
eq 8680 (atom 37887 = S*S) and eq 13985 (atom 39967, likewise).  Every other row of every model is
exactly affine, checked at steps 1, 2, 3 against all 49/53 knobs.  So the blast radius is those two
rows -- but they must be fixed everywhere they appear, which for the witness placement is inside the
139 OUTSIDE rows, where they were silently imposing the wrong constraint:
S = 0 at the witness, so a finite difference of S^2 measures (dS)^2, not dS.

Here: rebuild the witness model with both rows replaced by their linear cores, then re-run the
collateral kernel (widen) and the collateral-budget sweep (drop); and construct a real assignment
from the corrected |R|=13 model to check the OPT=6 the correction produced.
"""
import ast, json, itertools, os, re, sys, time
import ev, model
import optN
from optN import make, build, inner, WIT, POOL, atom_eqs
from widen import wide_knobs, int_kernel
import zsolve
from sqaudit import square_base, linear_core

HERE = os.path.dirname(os.path.abspath(__file__))
SQROWS = {}


def core(st, e):
    """linear core of row e: the inner form, stripped through any top-level square"""
    v, dep = linear_core(st, e)
    return v


def build_fixed(D):
    st = make(list(D))
    b0 = build(st)
    Rl = b0['R']
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows = list(Rl) + list(outside)
    b = [core(st, e) for e in rows]
    cols = []
    for Y in knobs:
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        cols.append([core(h, e) - b[i] for i, e in enumerate(rows)])
    k = len(knobs)
    M = [[cols[j][i] for j in range(k)] for i in range(len(rows))]
    # verify affineness of the fixed model
    bad = 0
    for j, Y in enumerate(knobs):
        h2 = st.clone().set_free({Y: st.fv.get(Y, 0) + 2})
        for i, e in enumerate(rows):
            if core(h2, e) - b[i] != 2 * M[i][j]:
                bad += 1
    return st, Rl, outside, rows, b, M, k, knobs, bad


def witness_analysis():
    print('=== witness placement, model rebuilt with the two square rows stripped ===', flush=True)
    st, Rl, outside, rows, b, M, k, knobs, bad = build_fixed(WIT)
    nR = len(Rl)
    print('  |R|=%d outside=%d knobs=%d ; non-affine (row,knob) pairs AFTER the fix: %d'
          % (nR, len(outside), k, bad), flush=True)
    nzout = [i for i in range(nR, len(rows)) if b[i] != 0]
    print('  outside rows nonzero at the witness after stripping: %d %s'
          % (len(nzout), [rows[i] for i in nzout]), flush=True)
    Mreg = [M[i] for i in range(nR)]
    breg = b[:nR]
    Mout = [M[i] for i in range(nR, len(rows))]
    opt_reg, rws, exh, _ = zsolve.max_zero_rows(Mreg, breg, k, nR)
    print('  region OPT with ALL %d knobs, no collateral limit: %d of %d (exh=%s)'
          % (k, opt_reg, nR, exh), flush=True)

    def g_of(keep):
        C = [Mout[i] for i in keep]
        K = int_kernel(C) if C else None
        if K is None:
            K = [[1 if a == bb else 0 for a in range(k)] for bb in range(k)]
        if not K:
            return 0, 0
        P = [[sum(Mreg[i][j] * v[j] for j in range(k)) for v in K] for i in range(nR)]
        o, _, _, _ = zsolve.max_zero_rows(P, breg, len(K), nR)
        return o, len(K)

    allout = list(range(len(outside)))
    g0, dim0 = g_of(allout)
    print('  |W|=0: kernel dim %d, g=%d, need %d -> %s' % (dim0, g0, nR - 6,
                                                           'BEATS' if g0 >= nR - 6 else 'no'),
          flush=True)
    best = nR - g0
    t0 = time.time()
    for dpt in (1, 2):
        need = (nR - 6) + dpt
        bg = -1
        bw = None
        for Wt in itertools.combinations(allout, dpt):
            keep = [i for i in allout if i not in set(Wt)]
            g, dim = g_of(keep)
            if g > bg:
                bg, bw = g, Wt
            if (nR - g) + dpt < best:
                best = (nR - g) + dpt
            if g >= need:
                print('  *** |W|=%d W=%s g=%d >= %d ***' % (dpt, list(Wt), g, need), flush=True)
        print('  |W|=%d: %d subsets exhaustive, max g=%d (W=%s), need %d -> %s ; best failing %d'
              ' (%.0fs)' % (dpt, len(list(itertools.combinations(allout, dpt))), bg,
                            list(bw) if bw else None, need,
                            'BEATS' if bg >= need else 'no', best, time.time() - t0), flush=True)
    print('  => corrected model, |W| <= 2: best failing = %d (score %d)' % (best, 39033 - best),
          flush=True)
    return best


def r13_realise():
    print('\n=== the |R|=13 correction: OPT 5 -> 6.  Which rows, and does it realise? ===',
          flush=True)
    for D in ([], [642], [29854], [642, 29854, 31864]):
        st = make(list(D))
        b0 = build(st)
        Rl, n, knobs = b0['R'], b0['n'], b0['knobs']
        b = [core(st, e) for e in Rl]
        M = []
        cols = []
        for Y in knobs:
            h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
            cols.append([core(h, e) - b[i] for i, e in enumerate(Rl)])
        M = [[cols[j][i] for j in range(n)] for i in range(len(Rl))]
        opt, rws, exh, _ = zsolve.max_zero_rows(M, b, n, len(Rl))
        eqs = [Rl[i] for i in rws]
        t = zsolve.witness_t(M, b, n, rws)
        ch = {}
        for j, Y in enumerate(knobs):
            if t and t[j]:
                ch[Y] = st.fv.get(Y, 0) + t[j]
        g = st.clone()
        if ch:
            g.set_free(ch)
        print('  D=%-22s OPT=%d rows=%s  8680 among them: %-5s  -> REAL score %d (model predicts %d)'
              % (str(list(D)), opt, eqs, 8680 in eqs, g.score(), 39033 - (len(Rl) - opt)),
              flush=True)


if __name__ == '__main__':
    r13_realise()
    witness_analysis()
