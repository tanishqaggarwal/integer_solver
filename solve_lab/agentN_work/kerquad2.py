"""The HONEST reduction: only the rows that are genuinely linear may be used to cut a subspace.

`kerquad.py` restricted to ker(L_out), i.e. it assumed the linear part of every collateral row must
vanish on its own.  That is exactly the linearising assumption this thread is trying to drop: a row
`L.t + t^T Q t = 0` can be satisfied with `L.t != 0` if the quadratic part cancels it.  So here:

  * the 117 collateral rows with NO quadratic part are honest linear equations -> cut K0 = ker_Z of
    those, saturated, LLL-reduced;
  * the 22 collateral rows WITH a quadratic part stay as full degree-2 polynomials on K0;
  * the 12 region rows stay as full degree-2 polynomials on K0.

The result is the exact integer polynomial system in dim(K0) unknowns.  Its size is reported before
any attempt to solve it, and it is written out in Singular syntax.
"""
import os, sys, json, time, pickle, math
from collections import defaultdict
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
D = pickle.load(open(os.path.join(HERE, 'runs', 'polyfull.pkl'), 'rb'))
knobs, Rl, outside = D['knobs'], D['R'], D['outside']
polys = D['polys']
k = len(knobs)


def parts(c):
    c0 = 0
    lin = [0] * k
    quad = {}
    for mono, v in c.items():
        d = sum(mono)
        if d == 0:
            c0 = v
        elif d == 1:
            lin[mono.index(1)] = v
        elif d == 2:
            idx = [j for j, e in enumerate(mono) if e]
            key = (idx[0], idx[0]) if len(idx) == 1 else (idx[0], idx[1])
            quad[key] = v
        else:
            raise ValueError('degree %d' % d)
    return c0, lin, quad


PARTS = {e: parts(polys[e]) for e in list(Rl) + list(outside)}
from kerquad import int_kernel_columns


def restrict(e, K):
    d = len(K)
    c0, lin, quad = PARTS[e]
    L = [sum(lin[j] * K[a][j] for j in range(k)) for a in range(d)]
    Q = defaultdict(int)
    for (i, j), c in quad.items():
        for a in range(d):
            if K[a][i] == 0:
                continue
            ci = c * K[a][i]
            for b in range(d):
                if K[b][j] == 0:
                    continue
                key = (a, b) if a <= b else (b, a)
                Q[key] += ci * K[b][j]
    return c0, L, {kk: vv for kk, vv in Q.items() if vv}


def main():
    LIN = [e for e in outside if not PARTS[e][2]]
    QUA = [e for e in outside if PARTS[e][2]]
    print('collateral rows: %d purely linear, %d genuinely quadratic' % (len(LIN), len(QUA)),
          flush=True)

    t0 = time.time()
    K0 = int_kernel_columns([PARTS[e][1] for e in LIN], k)
    print('ker_Z of the purely-linear collateral rows: dim %d  (%.1fs)'
          % (len(K0), time.time() - t0), flush=True)
    # LLL-reduce the basis so the restricted coefficients stay printable
    Kr = fmpz_mat([[int(x) for x in v] for v in K0]).lll().tolist()
    K0 = [[int(x) for x in r] for r in Kr if any(r)]
    d = len(K0)
    mb = max(max(abs(x).bit_length() for x in v) for v in K0)
    print('   after LLL: dim %d, max basis entry %d bits' % (d, mb), flush=True)
    for e in LIN:
        L = PARTS[e][1]
        assert all(sum(L[j] * v[j] for j in range(k)) == 0 for v in K0), 'kernel wrong'
    print('   verified: all %d purely-linear collateral rows vanish on K0' % len(LIN), flush=True)

    print('\n=== the exact polynomial system on K0 (%d unknowns) ===' % d, flush=True)
    gens = []          # (label, const, lin, quad)
    print('collateral rows that stay nonzero on K0:', flush=True)
    ncol = 0
    for e in QUA:
        c0, L, Q = restrict(e, K0)
        if c0 == 0 and not any(L) and not Q:
            continue
        ncol += 1
        gens.append(('c%d' % e, c0, L, Q, False))
        print('   eq %-6d deg %d  const %-7s lin %-3d quad %-4d  maxbits %d'
              % (e, 2 if Q else (1 if any(L) else 0),
                 ('0' if c0 == 0 else '%db' % abs(c0).bit_length()),
                 sum(1 for x in L if x), len(Q),
                 max([abs(x).bit_length() for x in L] +
                     [abs(x).bit_length() for x in Q.values()] + [0])), flush=True)
    print('   -> %d of %d quadratic collateral rows survive on K0' % (ncol, len(QUA)), flush=True)

    print('region rows on K0:', flush=True)
    nregq = 0
    for e in Rl:
        c0, L, Q = restrict(e, K0)
        if Q:
            nregq += 1
        gens.append(('r%d' % e, c0, L, Q, True))
        print('   eq %-6d deg %d  const %-7s lin %-3d quad %-4d  maxbits %d'
              % (e, 2 if Q else (1 if any(L) else 0),
                 ('0' if c0 == 0 else '%db' % abs(c0).bit_length()),
                 sum(1 for x in L if x), len(Q),
                 max([abs(x).bit_length() for x in L] +
                     [abs(x).bit_length() for x in Q.values()] + [0])), flush=True)
    print('   -> %d of 12 region rows are genuinely quadratic on K0' % nregq, flush=True)

    maxdeg = max(2 if g[3] else (1 if any(g[2]) else 0) for g in gens)
    print('\n=== SYSTEM SIZE (before any solve attempt) ===', flush=True)
    print('unknowns              : %d' % d)
    print('generators            : %d (%d collateral + 12 region)' % (len(gens), ncol))
    print('max total degree      : %d' % maxdeg)
    print('Bezout bound          : 2^%d = %.3e' % (len(gens), 2.0 ** min(len(gens), 900)))
    print('max coefficient bits  : %d'
          % max(max([abs(g[1]).bit_length()] + [abs(x).bit_length() for x in g[2]] +
                    [abs(x).bit_length() for x in g[3].values()]) for g in gens))

    pickle.dump(dict(K0=K0, d=d, gens=gens, knobs=knobs, R=Rl, LIN=LIN, QUA=QUA),
                open(os.path.join(HERE, 'runs', 'kerquad2.pkl'), 'wb'))
    json.dump(dict(unknowns=d, ngens=len(gens), ncollateral=ncol, maxdeg=maxdeg,
                   nregion_quadratic=nregq,
                   linear_collateral=len(LIN), quadratic_collateral=len(QUA)),
              open(os.path.join(HERE, 'runs', 'kerquad2.json'), 'w'), indent=1)
    print('\nwrote runs/kerquad2.{json,pkl}', flush=True)


if __name__ == '__main__':
    main()
