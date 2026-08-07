"""The witness region as an EXACT integer polynomial system, region rows AND collateral rows.

`polyexact.py` showed the 7 zero-collateral knobs give a system of total degree exactly 1 (the DAG
cone of those 7 knobs contains no downstream variable at all), while the 49 wide knobs give total
degree exactly 2.  So the wide model that `widen.py` / `drop.py` used — built by finite differences
with a step of 1 — is a *secant* of a genuinely quadratic map, not its linear part.

This script:
  1. builds the exact polynomials for the 12 region rows and every outside row the wide knobs touch;
  2. VERIFIES them against direct recomputation through the frame at random integer points;
  3. reports the exact monomial structure (which pairs of knobs multiply, coefficient sizes);
  4. separates the true linear part L from the step-1 secant S = L + diag(Q) that widen.py used;
  5. evaluates the exact collateral polynomials on the integer kernel of L, which is the only place
     a zero-collateral move can live, and measures whether the quadratic part destroys it.
"""
import os, sys, json, time, random, math
from collections import defaultdict

import ev, optN
from optN import make, build, WIT, POOL, fr, FREE, FR0, atom_eqs, _bits, inner
import frameB as FB
from polyexact import P

HERE = os.path.dirname(os.path.abspath(__file__))
eq_terms = ev.eq_terms


def wide_setup():
    st = make(WIT)
    b0 = build(st)
    Rl = b0['R']
    import widen
    knobs, outside = widen.wide_knobs(st, Rl, verbose=True)
    return st, Rl, knobs, outside


def exact_polys(st, rows, knobs):
    """Exact polynomial for the inner sum of every row in `rows`, in Z[t_0..t_{k-1}]."""
    P.NK = len(knobs)
    v = list(st.v)
    ns = {'v': v, '__builtins__': {}}
    aff, ck = set(), set()
    for j, Y in enumerate(knobs):
        v[Y] = P.var(j, st.fv.get(Y, 0))
        aff.update(fr.desc[Y])
        ck.update(fr.chk[Y])
    for u in sorted(aff, key=lambda u: fr.pos[u]):
        v[u] = eval(FB.DEFEXPR[u], ns)
    av = dict(st.av)
    for a in sorted(ck):
        av[a] = eval(FB.ACODE[a], ns)
    out = {}
    for e in rows:
        m, sq, tl = eq_terms[e]
        acc = P()
        for c, a in tl:
            x = av.get(a)
            if isinstance(x, P):
                acc = acc + x * c
            elif x:
                acc = acc + P.const(c * x)
        out[e] = acc
    return out


def evalP(pol, t):
    s = 0
    for k, c in pol.c.items():
        m = c
        for j, e in enumerate(k):
            if e:
                m *= t[j] ** e
        s += m
    return s


def verify(st, polys, rows, knobs, ntrial=6, seed=11):
    """Direct recomputation: set the knobs to base+t in the real frame and compare."""
    rnd = random.Random(seed)
    bad = 0
    for tr in range(ntrial):
        if tr == 0:
            t = [1] * len(knobs)
        elif tr == 1:
            t = [0] * len(knobs)
        else:
            t = [rnd.randint(-4, 4) for _ in knobs]
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + t[j] for j, Y in enumerate(knobs)})
        for e in rows:
            a = inner(h, e)
            b = evalP(polys[e], t)
            if a != b:
                bad += 1
                if bad <= 3:
                    print('   MISMATCH eq %d trial %d' % (e, tr), flush=True)
    print('verification: %d row-evaluations, %d mismatches' % (ntrial * len(rows), bad), flush=True)
    return bad == 0


def parts(pol, k):
    """(constant, linear vector, quadratic dict{(i,j):c}, higher dict) of an exact polynomial."""
    c0 = 0
    lin = [0] * k
    quad = {}
    hi = {}
    for mono, c in pol.c.items():
        d = sum(mono)
        if d == 0:
            c0 = c
        elif d == 1:
            lin[mono.index(1)] = c
        elif d == 2:
            idx = [j for j, e in enumerate(mono) if e]
            key = (idx[0], idx[0]) if len(idx) == 1 else (idx[0], idx[1])
            quad[key] = c
        else:
            hi[mono] = c
    return c0, lin, quad, hi


def main():
    st, Rl, knobs, outside = wide_setup()
    k = len(knobs)
    rows = list(Rl) + list(outside)
    print('region rows %d, outside rows %d, knobs %d' % (len(Rl), len(outside), k), flush=True)

    t0 = time.time()
    polys = exact_polys(st, rows, knobs)
    print('exact polynomials built in %.1fs' % (time.time() - t0), flush=True)

    ok = verify(st, polys, rows, knobs)
    if not ok:
        print('ABORT: polynomial model does not reproduce the frame'); return

    # --- structure --------------------------------------------------------------------------
    maxdeg = 0
    nq_rows = 0
    pairset = set()
    sq_vars = set()
    stats = []
    for e in rows:
        c0, lin, quad, hi = parts(polys[e], k)
        d = polys[e].deg()
        maxdeg = max(maxdeg, d)
        if quad or hi:
            nq_rows += 1
        for (i, j) in quad:
            if i == j:
                sq_vars.add(i)
            else:
                pairset.add((i, j))
        stats.append(dict(eq=e, region=e in set(Rl), deg=d, terms=polys[e].nterms(),
                          nquad=len(quad), nhigh=len(hi),
                          const_bits=(abs(c0).bit_length() if c0 else 0),
                          lin_bits=max([abs(x).bit_length() for x in lin] or [0]),
                          quad_bits=max([abs(c).bit_length() for c in quad.values()] or [0])))
    print('\n=== EXACT SYSTEM ===', flush=True)
    print('unknowns %d, equations %d (12 region + %d collateral), max total degree %d'
          % (k, len(rows), len(outside), maxdeg), flush=True)
    print('rows with a genuine degree>=2 part: %d of %d' % (nq_rows, len(rows)), flush=True)
    print('distinct quadratic monomials: %d cross-pairs + %d squares'
          % (len(pairset), len(sq_vars)), flush=True)
    qv = sorted(set([i for i, j in pairset] + [j for i, j in pairset]) | sq_vars)
    print('knobs appearing in ANY quadratic monomial: %d of %d -> %s'
          % (len(qv), k, [knobs[j] for j in qv]), flush=True)
    lin_only = [knobs[j] for j in range(k) if j not in set(qv)]
    print('knobs appearing only linearly: %d -> %s' % (len(lin_only), lin_only), flush=True)

    reg_stats = [s for s in stats if s['region']]
    out_stats = [s for s in stats if not s['region']]
    print('\nregion rows : %d of 12 quadratic; max coef bits const %d lin %d quad %d'
          % (sum(1 for s in reg_stats if s['nquad'] or s['nhigh']),
             max(s['const_bits'] for s in reg_stats), max(s['lin_bits'] for s in reg_stats),
             max(s['quad_bits'] for s in reg_stats)), flush=True)
    print('collateral  : %d of %d quadratic; max coef bits const %d lin %d quad %d'
          % (sum(1 for s in out_stats if s['nquad'] or s['nhigh']), len(out_stats),
             max(s['const_bits'] for s in out_stats), max(s['lin_bits'] for s in out_stats),
             max(s['quad_bits'] for s in out_stats)), flush=True)

    # --- the secant that widen.py used vs the true linear part ---------------------------------
    print('\n=== is widen.py\'s finite-difference matrix the true linear part? ===', flush=True)
    diff = 0
    for e in rows:
        c0, lin, quad, hi = parts(polys[e], k)
        for j in range(k):
            secant = lin[j] + quad.get((j, j), 0) + sum(c for m, c in hi.items()
                                                        if sum(m) == m[j])
            if secant != lin[j]:
                diff += 1
    print('(row,knob) entries where the step-1 secant differs from the true linear coefficient: %d'
          % diff, flush=True)

    json.dump(dict(knobs=knobs, R=Rl, outside=outside, maxdeg=maxdeg, nq_rows=nq_rows,
                   quad_pairs=len(pairset), quad_squares=len(sq_vars),
                   quad_knobs=[knobs[j] for j in qv], lin_only=lin_only,
                   secant_differs=diff, stats=stats),
              open(os.path.join(HERE, 'runs', 'polyfull.json'), 'w'), indent=1)

    # keep the polynomials for the solver stage
    import pickle
    pickle.dump(dict(knobs=knobs, R=Rl, outside=outside,
                     polys={e: polys[e].c for e in rows}),
                open(os.path.join(HERE, 'runs', 'polyfull.pkl'), 'wb'))
    print('\nwrote runs/polyfull.json and runs/polyfull.pkl', flush=True)


if __name__ == '__main__':
    main()
