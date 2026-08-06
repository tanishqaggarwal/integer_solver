"""IP #8 -- same integer program as IP #7 but with MODULAR PRE-SCREENING.

Exact integer solvability is expensive on 900-bit entries; solvability mod a random prime is a
cheap NECESSARY condition, so screen every candidate subset modulo two primes first and only
run the exact HNF solve on survivors.  This makes the subset search practical.
"""
import sys, os, json, itertools, time, collections, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip7 import atomval, load_raw, deltas
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
Q1, Q2 = (1 << 61) - 1, 2305843009213693907


def solvable_mod(M, rhs, q):
    m = len(M)
    n = len(M[0]) if m else 0
    A = [[M[i][j] % q for j in range(n)] + [rhs[i] % q] for i in range(m)]
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
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if A[i][n] and not any(A[i][:n]):
            return False
    return True


def build(v, verbose=True):
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
    core = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            core |= set(L.avars[a])
    dc, uc, collat = [], [], set()
    for u in sorted(core):
        d1 = deltas(v, AV, u, 1)
        if not d1:
            continue
        d2 = deltas(v, AV, u, 2)
        if not all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
            continue
        dc.append(d1)
        uc.append(u)
        collat |= {e for e in d1 if e not in set(FAIL)}
    comp = set()
    for e in collat:
        for a in L.eq_atoms[e][2]:
            comp |= set(L.avars[a])
    comp -= set(uc)
    dcomp, ucomp = [], []
    for u in sorted(comp):
        d1 = deltas(v, AV, u, 1)
        if not d1 or not (set(d1) & collat):
            continue
        d2 = deltas(v, AV, u, 2)
        if not all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
            continue
        if any(e not in collat and e not in set(FAIL) for e in d1):
            continue
        dcomp.append(d1)
        ucomp.append(u)
    used = uc + ucomp
    cols = dc + dcomp
    ROWS = list(FAIL) + sorted(collat)
    M = [[cols[j].get(e, 0) for j in range(len(cols))] for e in ROWS]
    rhs = [-eqs(e) for e in FAIL] + [0] * len(collat)
    if verbose:
        print(f"  failing={len(FAIL)} core={len(uc)} collateral={len(collat)} "
              f"compensators={len(ucomp)}  system {len(ROWS)}x{len(cols)}", flush=True)
    return v, FAIL, used, M, rhs, len(FAIL)


def search(v, maxallow=7, verbose=True):
    v, FAIL, used, M, rhs, nf = build(v, verbose)
    t0 = time.time()
    tested = screened = 0
    for allow in range(0, min(maxallow, nf) + 1):
        for combo in itertools.combinations(range(nf), allow):
            keep = [i for i in range(len(M)) if i not in set(combo)]
            Mk = [M[i] for i in keep]
            rk = [rhs[i] for i in keep]
            tested += 1
            if not solvable_mod(Mk, rk, Q1):
                continue
            if not solvable_mod(Mk, rk, Q2):
                continue
            screened += 1
            x = solve_int(Mk, rk)
            if x is None:
                continue
            snap = [v[u] for u in used]
            for j, u in enumerate(used):
                v[u] += x[j]
            AV2 = [atomval(a, v) for a in range(L.NA)]
            f2 = [e for e in range(L.NEQ)
                  if sum(c * AV2[a] for a, c in L.eq_atoms[e][2].items()) != 0]
            print(f"  allow={allow}: APPLIED -> failing={len(f2)} score={L.NEQ-len(f2)} "
                  f"({time.time()-t0:.0f}s, {tested} subsets tested, {screened} passed screening)",
                  flush=True)
            if len(f2) <= nf:
                return v, len(f2)
            for j, u in enumerate(used):
                v[u] = snap[j]
        if verbose:
            print(f"    allow={allow}: none ({tested} tested, {screened} passed, "
                  f"{time.time()-t0:.0f}s)", flush=True)
    return v, nf


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print(f"=== {os.path.basename(src)}  (IP with modular screening)")
    best = None
    prev = None
    for it in range(6):
        v, f = search(v)
        print(f"  it{it}: failing={f} score={L.NEQ-f}", flush=True)
        if best is None or f < best[0]:
            best = (f, [x for x in v])
        if f == 0 or f == prev:
            break
        prev = f
    print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
    if best[0] < 7:
        json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'ip8_best_named.json'), 'w'))
        print("saved data/ip8_best_named.json")
