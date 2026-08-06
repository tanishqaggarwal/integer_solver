"""THE RIGHT REFORMULATION -- minimum failing equations = SYNDROME DECODING.

Let K = {y : y^T M = 0 (mod p)} with basis B (rows), and w = B*rhs.
Dropping a row set D leaves the system solvable mod p  iff  every y in K vanishing on the
complement of D has y.rhs = 0, which is exactly

        w  in  span{ columns of B indexed by D }
   <=>  there is z supported on D with  B z = w.

So the MINIMUM NUMBER OF FAILING EQUATIONS (mod p, in this local system) is the minimum
weight of z with B z = w -- a syndrome-decoding problem over GF(p).

The checkpoint currently sits at weight 7.  Anything smaller beats it.
"""
import sys, os, json, time, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from ip8 import build
from dual import left_kernel_modp
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def rank_modp(cols, q):
    """rank of the list of column vectors"""
    if not cols:
        return 0
    n = len(cols[0])
    A = [list(c) for c in cols]
    r = 0
    for i in range(n):
        pr = None
        for j in range(r, len(A)):
            if A[j][i] % q:
                pr = j
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][i], -1, q)
        A[r] = [x * inv % q for x in A[r]]
        for j in range(len(A)):
            if j != r and A[j][i] % q:
                f = A[j][i]
                A[j] = [(A[j][k] - f * A[r][k]) % q for k in range(n)]
        r += 1
        if r == len(A):
            break
    return r


def in_span(cols, w, q):
    return rank_modp(cols + [w], q) == rank_modp(cols, q)


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print("===", os.path.basename(src))
    v, FAIL, used, M, rhs, nf = build(v)
    m = len(M)
    t0 = time.time()
    B, rank = left_kernel_modp(M, P)
    print(f"  rows={m} cols={len(M[0])} rank={rank} kernel dim={len(B)} ({time.time()-t0:.0f}s)")
    # columns of B, indexed by row of M
    cols = [[B[t][i] % P for t in range(len(B))] for i in range(m)]
    w = [sum(B[t][i] * rhs[i] for i in range(m)) % P for t in range(len(B))]
    print(f"  syndrome w nonzero: {any(w)}")
    nzc = [i for i in range(m) if any(cols[i])]
    print(f"  rows with a nonzero column (can contribute): {len(nzc)} of {m}")
    t0 = time.time()
    for wt in range(1, 5):
        found = None
        for D in itertools.combinations(nzc, wt):
            if in_span([cols[i] for i in D], w, P):
                found = D
                break
        if found:
            print(f"  MINIMUM WEIGHT = {wt}   D = {found}")
            print(f"    of those, currently-failing rows (idx < {nf}): {[i for i in found if i < nf]}")
            print(f"    => at most {wt} equations need fail in this local system "
                  f"(score up to {L.NEQ - wt})   ({time.time()-t0:.0f}s)")
            break
        print(f"    weight {wt}: none ({time.time()-t0:.0f}s)", flush=True)
