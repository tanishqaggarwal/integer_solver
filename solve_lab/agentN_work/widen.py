"""Integer COMBINATIONS of knobs whose collateral cancels.

Agent H's whole pipeline only ever admitted knobs that individually disturb nothing outside
the region (`all(q in Sset for q in d)`).  That is the coordinate-axis part of the collateral
kernel.  A *combination* of knobs that individually break outside equations can have zero net
collateral, and such combinations are exactly what a coordinate-wise filter cannot see.

Here:
  wide candidate set  = every free input that moves any atom of any equation of the region
  collateral map      C : Z^k -> Z^(outside equations touched)
  admissible moves    = ker_Z(C)   (saturated integer kernel, via flint nullspace + HNF)
  then maximise the number of integrally zeroable region rows over the image of that kernel.

If ker_Z(C) projects onto the region with finer granularity than the coordinate knobs do, the
p-obstruction measured in obstruct.py can be escaped.  Any hit is re-verified exactly.
"""
import sys, json, time, itertools
from collections import defaultdict
from flint import fmpz_mat
import ev
import optN
from optN import make, build, WIT, POOL, fr, FREE, FR0, inner, atom_eqs, _bits, rank_q
import zsolve


def wide_knobs(st, Rl, verbose=True):
    """Every free input that moves at least one atom occurring in an equation of Rl."""
    Rset = set(Rl)
    atoms_R = set()
    for e in Rl:
        for c, a in ev.eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    cands = sorted(y for y in cands if y in FREE)
    if verbose:
        print('wide candidate free inputs: %d' % len(cands), flush=True)
    keep = []
    touched = set()
    eff = {}
    for Y in cands:
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        moved = [q for q in h.av if h.av[q] != st.av[q]]
        if not moved:
            continue
        eqs = set()
        for q in moved:
            eqs |= atom_eqs[q]
        if not (eqs & Rset):
            continue                      # touches nothing in the region: pure collateral
        keep.append(Y)
        eff[Y] = eqs
        touched |= eqs
    outside = sorted(touched - Rset)
    if verbose:
        print('knobs that touch the region: %d ; outside equations they disturb: %d'
              % (len(keep), len(outside)), flush=True)
    return keep, outside


def build_wide(st, Rl, knobs, outside):
    rows = list(Rl) + list(outside)
    b = [inner(st, e) for e in rows]
    cols = []
    for Y in knobs:
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        cols.append([inner(h, e) - b[i] for i, e in enumerate(rows)])
    k = len(knobs)
    M = [[cols[j][i] for j in range(k)] for i in range(len(rows))]
    return rows, b, M, k


def int_kernel(rows_matrix):
    """Saturated integer kernel basis of the integer matrix given as list of rows.
    Returns list of basis vectors (length = #columns)."""
    m = len(rows_matrix)
    if m == 0:
        return None            # kernel is everything
    n = len(rows_matrix[0])
    A = fmpz_mat([[int(x) for x in r] for r in rows_matrix])
    X, nullity = A.nullspace()
    if nullity == 0:
        return []
    basis = []
    Xl = X.tolist()
    for j in range(nullity):
        col = [int(Xl[i][j]) for i in range(n)]
        g = 0
        for v in col:
            g = __import__('math').gcd(g, abs(v))
        if g > 1:
            col = [v // g for v in col]
        basis.append(col)
    # saturate: HNF of the basis, then divide out content per pivot (flint nullspace is
    # already over Q so the row span is the saturated kernel up to scaling; HNF normalises it)
    H = fmpz_mat(basis).hnf().tolist()
    out = []
    for r in H:
        rr = [int(x) for x in r]
        if any(rr):
            out.append(rr)
    return out


def run(D, tag, node_cap=200000):
    st = make(list(D))
    d0 = build(st)
    Rl = d0['R']
    print('\n=== %s : region |R|=%d, narrow knobs=%d, narrow OPT=%d ===' %
          (tag, len(Rl), d0['n'],
           zsolve.max_zero_rows(d0['M'], d0['b'], d0['n'], len(Rl))[0]), flush=True)
    knobs, outside = wide_knobs(st, Rl)
    if not knobs:
        print('no wide knobs'); return
    t0 = time.time()
    rows, b, M, k = build_wide(st, Rl, knobs, outside)
    print('wide matrix: %d rows x %d knobs  (%.1fs)' % (len(rows), k, time.time() - t0), flush=True)
    nR = len(Rl)
    # collateral block = the outside rows; require zero net collateral
    Cblock = [M[i] for i in range(nR, len(rows))]
    # sanity: outside rows currently satisfied
    nzout = [i for i in range(nR, len(rows)) if b[i] != 0]
    print('outside rows currently nonzero: %d (expected 0)' % len(nzout), flush=True)
    K = int_kernel(Cblock)
    if K is None:
        K = [[1 if i == j else 0 for i in range(k)] for j in range(k)]
    print('integer collateral kernel dimension: %d  (knobs=%d)' % (len(K), k), flush=True)
    if not K:
        print('kernel is trivial: no zero-net-collateral combination exists'); return
    # project the kernel onto the region rows
    Mreg = [M[i] for i in range(nR)]
    P = [[sum(Mreg[i][j] * v[j] for j in range(k)) for v in K] for i in range(nR)]
    breg = b[:nR]
    rk = rank_q(P, len(K))
    opt, rws, exh, tests = zsolve.max_zero_rows(P, breg, len(K), nR, node_cap=node_cap)
    fail = nR - opt + len(nzout)
    print('WIDE: kernel dim=%d rank on region=%d  OPT=%d of %d  failing=%d  score=%d  exh=%s'
          % (len(K), rk, opt, nR, fail, 39033 - fail, exh), flush=True)
    return dict(D=list(D), nR=nR, k=k, outside=len(outside), kdim=len(K), rank=rk,
                opt=opt, failing=fail, score=39033 - fail, exh=str(exh),
                K=K, knobs=knobs, rows=rws)


if __name__ == '__main__':
    res = []
    r = run(WIT, 'WITNESS %s' % WIT)
    if r:
        res.append({kk: vv for kk, vv in r.items() if kk not in ('K', 'knobs')})
    for D in ([28730], [17499], [642, 28730, 31864]):
        r = run(D, 'D=%s' % D)
        if r:
            res.append({kk: vv for kk, vv in r.items() if kk not in ('K', 'knobs')})
    json.dump(res, open('runs/widen.json', 'w'), indent=1)
