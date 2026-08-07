"""Reduce the exact quadratic system to the collateral kernel, and hand it to Singular.

The zero-collateral condition is 139 exact equations `L_r.t + t^T Q_r t = 0`, all with zero constant
term (the witness satisfies them).  Any solution must in particular satisfy the LINEAR part on the
subspace where the quadratic part vanishes; but exactly, the right first move is:

  ker_Z(L_out)  -- the integer lattice of moves with zero FIRST-ORDER collateral, computed with a
                   unimodular column reduction so it is genuinely saturated (every integer t with
                   L_out t = 0 is K s for an integer s), not merely a spanning set;

and then, on that lattice, the 139 collateral rows collapse to their pure quadratic parts, giving a
small system of quadratic forms in dim(K) variables.  That system is the exact obstruction the
linear model of widen.py/drop.py could not see.  Its solution set is what a wide-knob, zero-
collateral move actually has available.
"""
import os, sys, json, time, pickle, math
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE, 'runs', 'polyfull.pkl'), 'rb'))
knobs, Rl, outside = D['knobs'], D['R'], D['outside']
polys = D['polys']
k = len(knobs)
NK = k


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


def int_kernel_columns(A, n):
    """Saturated integer kernel of the integer matrix A (list of rows, n columns).

    Unimodular COLUMN reduction: keep U = I_n and apply the same column operations to U; the
    columns of U whose image column in A has become zero are a basis of ker_Z(A)."""
    m = len(A)
    cols = [[A[i][j] for i in range(m)] for j in range(n)]          # column-major
    U = [[1 if i == j else 0 for i in range(n)] for j in range(n)]  # U[j] = j-th column of U
    piv = 0
    for r in range(m):
        # find a column >= piv with a nonzero entry in row r
        idx = [j for j in range(piv, n) if cols[j][r] != 0]
        if not idx:
            continue
        # gcd-reduce those columns into one
        j0 = idx[0]
        cols[piv], cols[j0] = cols[j0], cols[piv]
        U[piv], U[j0] = U[j0], U[piv]
        for j in range(piv + 1, n):
            while cols[j][r] != 0:
                a, b = cols[piv][r], cols[j][r]
                if abs(a) > abs(b) or a == 0:
                    cols[piv], cols[j] = cols[j], cols[piv]
                    U[piv], U[j] = U[j], U[piv]
                    a, b = cols[piv][r], cols[j][r]
                q = b // a
                if q:
                    cols[j] = [x - q * y for x, y in zip(cols[j], cols[piv])]
                    U[j] = [x - q * y for x, y in zip(U[j], U[piv])]
                else:
                    break
        piv += 1
        if piv == n:
            break
    ker = [U[j] for j in range(n) if all(x == 0 for x in cols[j])]
    return ker


def main():
    print('knobs %d, region rows %d, collateral rows %d' % (k, len(Rl), len(outside)), flush=True)

    Lout = [PARTS[e][1] for e in outside]
    c_out = [PARTS[e][0] for e in outside]
    assert all(c == 0 for c in c_out), 'collateral rows must be satisfied at the witness'

    t0 = time.time()
    K = int_kernel_columns(Lout, k)
    print('saturated integer kernel of the TRUE linear collateral part: dim %d  (%.1fs)'
          % (len(K), time.time() - t0), flush=True)
    # exactness check
    bad = 0
    for r, row in enumerate(Lout):
        for v in K:
            if sum(row[j] * v[j] for j in range(k)) != 0:
                bad += 1
    print('   L_out . K == 0 check: %d violations' % bad, flush=True)

    # widen.py's secant matrix, for the blast-radius comparison
    Sout = []
    for e in outside:
        c0, lin, quad = PARTS[e]
        Sout.append([lin[j] + quad.get((j, j), 0) for j in range(k)])
    Ks = int_kernel_columns(Sout, k)
    print('kernel of widen.py\'s STEP-1 SECANT matrix: dim %d  (widen.py reported 14)'
          % len(Ks), flush=True)

    # do the 7 narrow (zero-collateral) knobs sit in the kernel?
    import optN
    st = optN.make(optN.WIT)
    narrow = optN.build(st)['knobs']
    ni = [knobs.index(Y) for Y in narrow if Y in knobs]
    print('narrow knobs present among the 49: %d of %d' % (len(ni), len(narrow)), flush=True)
    for j in ni:
        e = [0] * k
        e[j] = 1
        assert all(sum(row[jj] * e[jj] for jj in range(k)) == 0 for row in Lout), \
            'narrow knob %d has nonzero true linear collateral' % knobs[j]
    print('   all narrow knobs verified to lie in ker(L_out) exactly', flush=True)

    # --- restrict every row to the kernel ------------------------------------------------------
    d = len(K)

    def restrict(e):
        c0, lin, quad = PARTS[e]
        # t = sum_a s_a K[a] ; linear part
        L = [sum(lin[j] * K[a][j] for j in range(k)) for a in range(d)]
        Q = defaultdict(int)
        for (i, j), c in quad.items():
            for a in range(d):
                if K[a][i] == 0:
                    continue
                for b in range(d):
                    if K[b][j] == 0:
                        continue
                    key = (a, b) if a <= b else (b, a)
                    Q[key] += c * K[a][i] * K[b][j]
        Q = {kk: vv for kk, vv in Q.items() if vv}
        return c0, L, Q

    print('\n=== the 139 collateral rows restricted to ker(L_out) ===', flush=True)
    coll = []
    nonzero_quadratic = 0
    for e in outside:
        c0, L, Q = restrict(e)
        assert c0 == 0 and all(x == 0 for x in L), 'linear part did not vanish on the kernel'
        if Q:
            nonzero_quadratic += 1
            coll.append((e, Q))
    print('rows that vanish identically on the kernel: %d' % (len(outside) - nonzero_quadratic),
          flush=True)
    print('rows that survive as NONZERO QUADRATIC FORMS in %d variables: %d'
          % (d, nonzero_quadratic), flush=True)
    if coll:
        print('   their monomial counts: %s' % sorted(len(Q) for e, Q in coll), flush=True)
        print('   max coefficient bits: %d'
              % max(max(abs(v).bit_length() for v in Q.values()) for e, Q in coll), flush=True)

    print('\n=== the 12 region rows restricted to ker(L_out) ===', flush=True)
    reg = []
    for e in Rl:
        c0, L, Q = restrict(e)
        reg.append((e, c0, L, Q))
        print('   eq %-6d const %-6s linterms %-3d quadterms %-4d'
              % (e, ('0' if c0 == 0 else '%db' % abs(c0).bit_length()),
                 sum(1 for x in L if x), len(Q)), flush=True)

    pickle.dump(dict(K=K, d=d, knobs=knobs, narrow=narrow, coll=coll, reg=reg,
                     outside=outside, R=Rl),
                open(os.path.join(HERE, 'runs', 'kerquad.pkl'), 'wb'))
    json.dump(dict(kernel_dim=d, secant_kernel_dim=len(Ks),
                   collateral_quadric_count=nonzero_quadratic,
                   collateral_quadric_terms=sorted(len(Q) for e, Q in coll),
                   region=[dict(eq=e, const_bits=(abs(c0).bit_length() if c0 else 0),
                                lin=sum(1 for x in L if x), quad=len(Q))
                           for e, c0, L, Q in reg]),
              open(os.path.join(HERE, 'runs', 'kerquad.json'), 'w'), indent=1)
    print('\nwrote runs/kerquad.{json,pkl}', flush=True)


if __name__ == '__main__':
    main()
