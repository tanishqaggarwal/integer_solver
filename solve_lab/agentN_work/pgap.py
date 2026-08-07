"""Price rk_p(M) against rk_Q(M) across configurations: is the mod-p gap an invariant?

The barrier at the deliverable is a mod-p rank deficiency of the region response matrix on the
zero-collateral lattice: rk_Q(M) = 7 but rk_p(M) = 3, and rk_p([M|b]) = 4, so the system is
inconsistent mod p.  `pgrow.py` showed the gap `rk_p([M|b]) - rk_p(M)` stays exactly 1 as the
lattice is enlarged by paying collateral.  This asks the sharper question: does the gap depend on
the CONFIGURATION?

Configurations priced:
  A. the 16 detach states, which are the whole 2^65 detach lattice by proof;
  B. every independently verified assignment on disk, loaded as a base state.

For each: the complete knob set (exact syntactic support, not a step-1 filter), the exact
saturation loop, then rk_Q / rk_p of M and [M|b], the gap, and the exhaustive integer OPT.
One rank computation per configuration — no sweeps.

Usage: python3 pgap.py
"""
import os, sys, json, time, glob, itertools
from collections import defaultdict
from fractions import Fraction
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN, zsolve
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits
from frameB import State
from polyexact import P
from polyfull import exact_polys
from kerquad import int_kernel_columns
from sqaudit import square_base
import frameB as FB
import re as _re
_VR = _re.compile(r'x_(\d+)')
_BASECODE = {}


_MISS = object()


def base_code(a):
    c = _BASECODE.get(a, _MISS)
    if c is _MISS:
        sb = square_base(a)
        _BASECODE[a] = c = (compile(_VR.sub(r'v[\1]', sb), '<sb>', 'eval') if sb else None)
    return c

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663
eq_terms = ev.eq_terms


def lll(K):
    if not K:
        return []
    R = fmpz_mat([[int(x) for x in v] for v in K]).lll().tolist()
    return [[int(x) for x in r] for r in R if any(r)]


def rank_mod(rows, ncol, p):
    A = [[x % p for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(ncol)]
        r += 1
    return r


def rank_q(rows, ncol):
    A = [[Fraction(x) for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(ncol)]
        r += 1
    return r


def price(st, tag):
    """Full exact pricing of one configuration.  Returns a dict or None if the region is empty."""
    NZ = set(st.nz())
    if not NZ:
        return dict(tag=tag, note='no nonzero atoms — already perfect?')
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    Rl = sorted(R)
    outside_fail = len([e for e in st.fails if e not in R])

    atoms_R = set()
    for e in Rl:
        for c, a in eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    cands = sorted(y for y in cands if y in FREE)
    k = len(cands)
    touched = set()
    for Y in cands:
        for a in fr.chk.get(Y, []):
            touched |= atom_eqs[a]
    outside = sorted(touched - R)

    # symbolic environment: knobs symbolic, everything else the base state's integer value
    P.NK = k
    vv = list(st.v)
    ns = {'v': vv, '__builtins__': {}}
    aff, ck = set(), set()
    for j, Y in enumerate(cands):
        vv[Y] = P.var(j, st.fv.get(Y, 0))
        aff.update(fr.desc[Y])
        ck.update(fr.chk[Y])
    for u in sorted(aff, key=lambda u: fr.pos[u]):
        vv[u] = eval(FB.DEFEXPR[u], ns)
    av = dict(st.av)
    for a in sorted(ck):
        av[a] = eval(FB.ACODE[a], ns)

    def rowpoly(e):
        """exact polynomial whose vanishing is equivalent to equation e being satisfied.
        A row that is a single top-level SQUARE atom is replaced by that atom's BASE: the
        equation is a power of the base, so `row = 0` iff `base = 0`, and truncating the
        square instead (as a linear model does) is exactly the error T caught in eq 8680."""
        m, sq, tl = eq_terms[e]
        live_tl = [(c, a) for c, a in tl if isinstance(av.get(a), P) or av.get(a)]
        if len(live_tl) == 1 and base_code(live_tl[0][1]) is not None:
            x = eval(base_code(live_tl[0][1]), ns)
            return x if isinstance(x, P) else P.const(x), True
        acc = P()
        for c, a in tl:
            x = av.get(a)
            if isinstance(x, P):
                acc = acc + x * c
            elif x:
                acc = acc + P.const(c * x)
        return acc, False

    polys = {}
    rooted = []
    for e in Rl + outside:
        pol, wasq = rowpoly(e)
        polys[e] = pol
        if wasq:
            rooted.append(e)
    live = [e for e in outside if polys[e].c]
    # collateral rows must be satisfied at the base state
    broken = [e for e in live if polys[e].c.get((0,) * k, 0) != 0]

    def restrict(pol, K):
        dd = len(K)
        P.NK = dd
        T = []
        for j in range(k):
            c = {}
            for a in range(dd):
                if K[a][j]:
                    m = [0] * dd
                    m[a] = 1
                    c[tuple(m)] = K[a][j]
            T.append(P(c))
        acc = P()
        for mono, cf in pol.c.items():
            term = P.const(cf)
            for j, ex in enumerate(mono):
                for _ in range(ex):
                    term = term * T[j]
            acc = acc + term
        return acc

    K = [[1 if i == j else 0 for i in range(k)] for j in range(k)]
    live_o = [e for e in live if e not in set(broken)]
    maxdeg = max([polys[e].deg() for e in Rl + live] or [0])
    for _ in range(14):
        cur = {e: restrict(polys[e], K) for e in live_o}
        linrows = [e for e in live_o if cur[e].deg() == 1]
        nl = [e for e in live_o if cur[e].deg() >= 2]
        if not linrows:
            live_o = nl
            break
        A = []
        for e in linrows:
            v = [0] * len(K)
            for mono, c in cur[e].c.items():
                v[mono.index(1)] = c
            A.append(v)
        Kn = int_kernel_columns(A, len(K))
        if not Kn:
            return dict(tag=tag, note='lattice collapsed to 0')
        K = lll([[sum(u[a] * K[a][j] for a in range(len(u))) for j in range(k)]
                 for u in lll(Kn)])
        live_o = nl
    dead = set()
    resid = []
    for e in live_o:
        c = restrict(polys[e], K).c
        if not c:
            continue
        if len(c) == 1:
            mono = list(c)[0]
            idx = [j for j, x in enumerate(mono) if x]
            if len(idx) == 1 and mono[idx[0]] == 2:
                dead.add(idx[0])
                continue
        resid.append(e)
    free = [a for a in range(len(K)) if a not in dead]

    M, b, nonlin = [], [], 0
    for e in Rl:
        c = restrict(polys[e], K).c
        if any(sum(m) > 1 for m in c):
            nonlin += 1
        row = [0] * len(K)
        c0 = 0
        for mono, v in c.items():
            if sum(mono) == 0:
                c0 = v
            elif sum(mono) == 1:
                row[mono.index(1)] = v
        M.append([row[a] for a in free])
        b.append(c0)
    n = len(free)
    aug = [M[i] + [b[i]] for i in range(len(M))]
    rQ, rQa = rank_q(M, n), rank_q(aug, n + 1)
    rP, rPa = rank_mod(M, n, Pp), rank_mod(aug, n + 1, Pp)
    opt, rws, exh, nd = zsolve.max_zero_rows(M, b, n, len(M), node_cap=2000000)
    fail = len(Rl) - opt + outside_fail
    return dict(tag=tag, R=len(Rl), knobs=k, lattice=n, maxdeg=maxdeg, rooted=rooted,
                rk_Q=rQ, rk_Q_aug=rQa, gap_Q=rQa - rQ,
                rk_p=rP, rk_p_aug=rPa, gap_p=rPa - rP,
                opt=opt, exh=bool(exh), outside=outside_fail,
                failing=fail, score=39033 - fail,
                residual_nonlinear=len(resid), region_nonlinear=nonlin,
                broken_collateral=len(broken))


HDR = ('%-34s %-4s %-6s %-8s %-6s %-7s %-6s %-7s %-6s %-5s %-6s'
       % ('configuration', '|R|', 'knobs', 'lattice', 'rk_Q', 'gap_Q', 'rk_p', 'gap_p', 'OPT',
          'fail', 'score'))


def show(r):
    if 'rk_Q' not in r:
        print('%-34s %s' % (r['tag'], r.get('note')), flush=True)
        return
    print('%-34s %-4d %-6d %-8d %-6d %-7d %-6d %-7d %-6d %-5d %-6d'
          % (r['tag'], r['R'], r['knobs'], r['lattice'], r['rk_Q'], r['gap_Q'],
             r['rk_p'], r['gap_p'], r['opt'], r['failing'], r['score']), flush=True)


def main():
    out = []
    print('=== A. the 16 detach states (the whole 2^65 lattice, by proof) ===', flush=True)
    print(HDR, flush=True)
    W4 = [642, 28730, 29854, 31864]
    for m in range(16):
        D = [W4[i] for i in range(4) if m >> i & 1]
        t0 = time.time()
        r = price(make(D), 'D=%s' % (D if D else '[]'))
        r['secs'] = round(time.time() - t0, 1)
        show(r)
        out.append(r)

    print('\n=== B. independently verified assignments on disk, as base configurations ===',
          flush=True)
    print(HDR, flush=True)
    files = sorted(glob.glob(os.path.join(HERE, '*.json')) +
                   glob.glob(os.path.join(HERE, '..', 'best', '*.json')))
    for f in files:
        try:
            W = json.load(open(f))
        except Exception:
            continue
        if not isinstance(W, dict) or len(W) < 1000:
            continue
        try:
            v = [0] * 38748
            for kk, val in W.items():
                v[int(kk[2:]) if str(kk).startswith('x_') else int(kk)] = int(val)
            fv = {u: v[u] for u in fr.free if v[u] != 0}
            st = State(fr, fv)
        except Exception as e:
            print('%-34s load failed: %s' % (os.path.basename(f), e), flush=True)
            continue
        tag = os.path.basename(f)[:32]
        t0 = time.time()
        try:
            r = price(st, tag)
        except Exception as e:
            print('%-34s price failed: %s' % (tag, e), flush=True)
            continue
        r['secs'] = round(time.time() - t0, 1)
        r['file'] = f
        r['frame_score'] = st.score()
        show(r)
        out.append(r)

    gaps = sorted(set(r['gap_p'] for r in out if 'gap_p' in r))
    gq = sorted(set(r['gap_Q'] for r in out if 'gap_Q' in r))
    print('\n=== VERDICT ===', flush=True)
    print('configurations priced        : %d' % len([r for r in out if 'gap_p' in r]), flush=True)
    print('distinct mod-p gaps observed : %s' % gaps, flush=True)
    print('distinct over-Q gaps observed: %s  (0 = solvable over Q)' % gq, flush=True)
    best = max((r for r in out if 'score' in r), key=lambda r: r['score'], default=None)
    if best:
        print('best score over all priced configurations: %d (%s)'
              % (best['score'], best['tag']), flush=True)
    json.dump(out, open(os.path.join(HERE, 'runs', 'pgap.json'), 'w'), indent=1)
    print('wrote runs/pgap.json', flush=True)


if __name__ == '__main__':
    main()
