"""Build the CLOSED region instead of the artificially filtered one.

ip8.build only admitted compensators disturbing nothing outside the region -- that threw away
126 of 176 candidates.  Here we instead let a compensator in and ADD the equations it disturbs
as new rows, iterating to a fixpoint (with a size cap).  Bigger system, but a far bigger column
space, which is what the obstruction is measured against.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip7 import atomval, load_raw
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def build_closed(v, maxrows=900, verbose=True):
    AV = [atomval(a, v) for a in range(L.NA)]

    def eqs(e):
        return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())

    def delta(u, step):
        old = v[u]
        v[u] = old + step
        d = collections.defaultdict(int)
        for a in L.var_atoms[u]:
            dv = atomval(a, v) - AV[a]
            if dv:
                for e, c in L.atom2eq.get(a, {}).items():
                    d[e] += c * dv
        v[u] = old
        return d

    FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
    ROWS = set(FAIL)
    VARS = {}                       # var -> delta dict
    frontier = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            frontier |= set(L.avars[a])
    seen = set()
    for it in range(12):
        added = 0
        for u in sorted(frontier):
            if u in seen:
                continue
            seen.add(u)
            d1 = delta(u, 1)
            if not d1:
                continue
            d2 = delta(u, 2)
            if not all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
                continue          # not exact-linear
            if len(ROWS | set(d1)) > maxrows:
                continue          # would blow the cap
            VARS[u] = d1
            ROWS |= set(d1)
            added += 1
        if verbose:
            print(f"    closure it{it}: rows={len(ROWS)} vars={len(VARS)} (+{added})", flush=True)
        if not added:
            break
        frontier = set()
        for e in ROWS:
            for a in L.eq_atoms[e][2]:
                frontier |= set(L.avars[a])
        frontier -= seen
    ROWS = sorted(ROWS)
    used = sorted(VARS)
    M = [[VARS[u].get(e, 0) for u in used] for e in ROWS]
    rhs = [-eqs(e) for e in ROWS]
    return v, FAIL, ROWS, used, M, rhs


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 900
    v = load_raw(src)
    print(f"=== {os.path.basename(src)}  (closed region, cap {cap})")
    t0 = time.time()
    v, FAIL, ROWS, used, M, rhs = build_closed(v, maxrows=cap)
    print(f"  FAIL={len(FAIL)}  closed system {len(ROWS)} x {len(used)}  ({time.time()-t0:.0f}s)",
          flush=True)
    t1 = time.time()
    x = solve_int(M, rhs)
    print(f"  full system integer-solvable: {'YES' if x else 'no'}  ({time.time()-t1:.0f}s)")
    if x:
        for j, u in enumerate(used):
            v[u] += x[j]
        AV = [atomval(a, v) for a in range(L.NA)]
        F = [e for e in range(L.NEQ)
             if sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) != 0]
        print(f"  APPLIED -> failing={len(F)} score={L.NEQ-len(F)}")
        json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'closed_solved.json'), 'w'))
        print("  saved data/closed_solved.json")
