"""IP #14 -- can the failing RIGHT-HAND SIDE be made p-divisible?

The invariant always contains exactly one factor of p.  A necessary (and, given the small
cofactors are CRT-clearable, essentially sufficient) step to a full solve is to reach a state
where the failing equation values are already = 0 (mod p).

That question is LINEAR OVER GF(p) -- fast, unlike the integer HNF work:

        find k with   b + G k  =  0   (mod p)   on the failing equations
                      G k      =  0   (mod p)   on the satisfied ones
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import atomval, load_raw, deltas
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def gf_solve(M, rhs, q):
    """solve M x = rhs over GF(q); return x or None"""
    m = len(M)
    n = len(M[0]) if m else 0
    A = [[M[i][j] % q for j in range(n)] + [rhs[i] % q] for i in range(m)]
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, m):
            if A[i][c]:
                pr = i
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], -1, q)
        A[r] = [x * inv % q for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][k] - f * A[r][k]) % q for k in range(n + 1)]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if A[i][n] % q and not any(A[i][:n]):
            return None
    x = [0] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x


def analyse(path, name):
    v = load_raw(path)
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
    print(f"=== {name}: failing={len(FAIL)}")
    print(f"    failing values = 0 mod p : {sum(1 for e in FAIL if eqs(e) % P == 0)} of {len(FAIL)}")
    core = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            core |= set(L.avars[a])
    # mod p we can use EVERY variable that moves the region, incl. the quadratic ones
    E = set(FAIL)
    for u in core:
        for a in L.var_atoms[u]:
            E |= set(L.atom2eq.get(a, {}))
    cand = set()
    for e in E:
        for a in L.eq_atoms[e][2]:
            cand |= set(L.avars[a])
    cand = sorted(cand)
    E = sorted(E)
    idx = {e: i for i, e in enumerate(E)}
    cols, used = [], []
    for u in cand:
        d = deltas(v, AV, u, 1)
        col = [d.get(e, 0) % P for e in E]
        if any(col):
            cols.append(col)
            used.append(u)
    M = [[cols[j][i] for j in range(len(used))] for i in range(len(E))]
    rhs = [(-eqs(e)) % P for e in E]
    print(f"    region {len(E)} equations, {len(used)} variables -> GF(p) system", flush=True)
    t0 = time.time()
    x = gf_solve(M, rhs, P)
    print(f"    solvable over GF(p): {x is not None}  ({time.time()-t0:.0f}s)")
    if x is None:
        # how many failing equations CAN be made p-divisible?
        FS = set(FAIL)
        keep = [i for i, e in enumerate(E) if e not in FS]
        x2 = gf_solve([M[i] for i in keep], [rhs[i] for i in keep], P)
        print(f"    (keeping only the satisfied equations at 0 mod p: "
              f"{'solvable' if x2 is not None else 'not solvable'})")
    return x is not None


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    for rel, nm in [(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'), 'checkpoint 39026'),
                    (os.path.join(HERE, 'data', 'finish3_named.json'), 's11 best 39018'),
                    (os.path.join(HERE, 'data', 'closehit2.json'), 'closehit2 39005')]:
        try:
            analyse(rel, nm)
        except Exception as e:
            print(nm, 'error', e)
