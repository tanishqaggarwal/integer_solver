"""IP #15 -- two-stage solve suggested by the invariant.

Stage A (mod p, fast):   drive EVERY equation in the region to 0 mod p.  The GF(p) system is
                         solvable (IP #14), and the map is polynomial, so iterate.
Stage B (over Z):        with the right-hand side now p-divisible, the single factor of p in
                         the invariant is absorbed -- run the exact integer solve.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip7 import atomval, load_raw, deltas
from ip14 import gf_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def region_of(v, AV, FAIL):
    core = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            core |= set(L.avars[a])
    E = set(FAIL)
    for u in core:
        for a in L.var_atoms[u]:
            E |= set(L.atom2eq.get(a, {}))
    cand = set()
    for e in E:
        for a in L.eq_atoms[e][2]:
            cand |= set(L.avars[a])
    return sorted(E), sorted(cand)


def stage_a(v, rounds=8, verbose=True):
    """drive every region equation to 0 mod p"""
    for rnd in range(rounds):
        AV = [atomval(a, v) for a in range(L.NA)]
        def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
        FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
        badp = [e for e in range(L.NEQ) if eqs(e) % P != 0]
        if verbose:
            print(f"  A{rnd}: failing={len(FAIL)}  nonzero mod p={len(badp)}", flush=True)
        if not badp:
            return v, True
        E, cand = region_of(v, AV, badp)
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
        t0 = time.time()
        x = gf_solve(M, rhs, P)
        if x is None:
            if verbose:
                print("     GF(p) system unsolvable at this point")
            return v, False
        for j, u in enumerate(used):
            v[u] = (v[u] + x[j])
        if verbose:
            print(f"     applied GF(p) step over {len(used)} variables ({time.time()-t0:.0f}s)",
                  flush=True)
    AV = [atomval(a, v) for a in range(L.NA)]
    badp = [e for e in range(L.NEQ)
            if sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) % P != 0]
    return v, not badp


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'data', 'closehit2.json')
    v = load_raw(src)
    AV = [atomval(a, v) for a in range(L.NA)]
    F0 = [e for e in range(L.NEQ) if sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) != 0]
    print(f"=== {os.path.basename(src)}: start failing={len(F0)} score={L.NEQ-len(F0)}")
    t0 = time.time()
    v, ok = stage_a(v)
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    F = [e for e in range(L.NEQ) if eqs(e) != 0]
    badp = [e for e in range(L.NEQ) if eqs(e) % P != 0]
    print(f"STAGE A done ({time.time()-t0:.0f}s): all-region 0 mod p = {ok}")
    print(f"  failing={len(F)} score={L.NEQ-len(F)}  nonzero mod p={len(badp)}")
    if ok and F:
        print("  every failing value is now a multiple of p -> exact integer solve")
        E, cand = region_of(v, AV, F)
        idx = {e: i for i, e in enumerate(E)}
        cols, used = [], []
        for u in cand:
            d1 = deltas(v, AV, u, 1)
            if not d1:
                continue
            d2 = deltas(v, AV, u, 2)
            if all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
                cols.append(d1)
                used.append(u)
        M = [[cols[j].get(e, 0) for j in range(len(cols))] for e in E]
        rhs = [-eqs(e) for e in E]
        print(f"  system {len(E)} x {len(cols)}", flush=True)
        t1 = time.time()
        x = solve_int(M, rhs)
        print(f"  exact integer solve: {'FOUND' if x else 'none'} ({time.time()-t1:.0f}s)")
        if x:
            for j, u in enumerate(used):
                v[u] += x[j]
            AV2 = [atomval(a, v) for a in range(L.NA)]
            F2 = [e for e in range(L.NEQ)
                  if sum(c * AV2[a] for a, c in L.eq_atoms[e][2].items()) != 0]
            print(f"  FINAL failing={len(F2)} score={L.NEQ-len(F2)}")
            json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                      open(os.path.join(HERE, 'data', 'ip15_named.json'), 'w'))
            print("  saved data/ip15_named.json")
