"""Zero the whole region and pay the collateral: is the trade better than 7?

With the full wide knob set (no zero-collateral filter) every region row is integrally
zeroable.  The solution set is t0 + ker_Z(M_region).  Search that coset for the member with
the least collateral damage, then evaluate the REAL score by exact re-evaluation (the linear
model is only a model).
"""
import json, time, itertools, random, sys
from flint import fmpz_mat
import zsolve
import optN
from optN import make, build, WIT, rank_q, inner
from widen import wide_knobs, build_wide, int_kernel


def apply_t(st, knobs, t):
    ch = {}
    for j, Y in enumerate(knobs):
        if t[j]:
            ch[Y] = st.fv.get(Y, 0) + t[j]
    g = st.clone()
    if ch:
        g.set_free(ch)
    return g


def analyse(D, tag, tries=4000, seed=1):
    random.seed(seed)
    st = make(list(D))
    d0 = build(st)
    Rl = d0['R']
    nR = len(Rl)
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows, b, M, k = build_wide(st, Rl, knobs, outside)
    Mreg = [M[i] for i in range(nR)]
    breg = b[:nR]
    print('\n=== %s ===  |R|=%d wide knobs=%d outside touched=%d  base score=%d'
          % (tag, nR, k, len(outside), st.score()), flush=True)
    t0 = zsolve.witness_t(Mreg, breg, k, range(nR))
    if t0 is None:
        print('  region not fully zeroable'); return None
    print('  t0 found; max |t| bit-length = %d, nonzero knobs = %d'
          % (max((x.bit_length() for x in t0), default=0), sum(1 for x in t0 if x)), flush=True)
    # predicted collateral under the linear model
    def model_damage(t):
        dmg = 0
        for i in range(nR, len(rows)):
            v = b[i] + sum(t[j] * M[i][j] for j in range(k))
            if v:
                dmg += 1
        return dmg
    print('  linear-model collateral of t0: %d of %d outside rows  => model failing=%d'
          % (model_damage(t0), len(outside), model_damage(t0)), flush=True)
    # the real thing
    g = apply_t(st, knobs, t0)
    print('  ACTUAL score after applying t0: %d  (failing %d)' % (g.score(), len(g.fails)), flush=True)
    best = (g.score(), t0, g)
    # search the coset t0 + ker(Mreg) for less damage
    K = int_kernel(Mreg)
    print('  ker_Z(M_region) dimension: %d' % (len(K) if K is not None else k), flush=True)
    if K:
        # greedy/random coset search on the linear model, then evaluate the good ones
        cur = list(t0)
        curd = model_damage(cur)
        print('  coset search over %d kernel directions ...' % len(K), flush=True)
        improved = True
        rounds = 0
        while improved and rounds < 30:
            improved = False
            rounds += 1
            for v in K:
                for s in (1, -1):
                    cand = [cur[j] + s * v[j] for j in range(k)]
                    dmg = model_damage(cand)
                    if dmg < curd:
                        cur, curd = cand, dmg
                        improved = True
        print('  best model collateral after coset descent: %d' % curd, flush=True)
        g2 = apply_t(st, knobs, cur)
        print('  ACTUAL score for the descended point: %d' % g2.score(), flush=True)
        if g2.score() > best[0]:
            best = (g2.score(), cur, g2)
    return best


if __name__ == '__main__':
    out = []
    for D, tag in ((WIT, 'WITNESS %s' % WIT), ([28730], 'D=[28730]'), ([17499], 'D=[17499]')):
        r = analyse(D, tag)
        if r and r[0] > 39026:
            sc, t, g = r
            path = 'N_%d.json' % sc
            json.dump({('x_%d' % i): g.v[i] for i in range(38748) if g.v[i] != 0},
                      open(path, 'w'))
            print('*** WROTE %s (score %d) ***' % (path, sc), flush=True)
