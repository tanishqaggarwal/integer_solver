"""IP #17 -- the two-stage solve done correctly.

IP #15 drove every region equation to 0 mod p but that turned satisfied equations into nonzero
multiples of p (28 -> 6097 failing).  The fix: move only inside the INTEGER KERNEL of the
satisfied equations, so they stay EXACTLY zero, and inside that kernel drive the failing values
to 0 mod p.  Then the single factor of p in the invariant is absorbed and the exact integer
solve should close.

  Stage A : k in ker_Z(G_keep)   with   b_fail + G_fail k = 0 (mod p)     [GF(p) inside the kernel]
  Stage B : exact integer solve of the now p-divisible residual
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip3 import int_kernel
from ip7 import atomval, load_raw, deltas
from ip14 import gf_solve
from ip15 import region_of
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def build_linear(v, AV, FAIL):
    E, cand = region_of(v, AV, FAIL)
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
    return E, idx, cols, used


def run(path, name, verbose=True):
    v = load_raw(path)
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
    print(f"=== {name}: failing={len(FAIL)} score={L.NEQ-len(FAIL)}")
    t0 = time.time()
    E, idx, cols, used = build_linear(v, AV, FAIL)
    FS = set(FAIL)
    KEEP = [e for e in E if e not in FS]
    print(f"  region {len(E)} eqs, {len(used)} exact-linear vars, keep {len(KEEP)} "
          f"({time.time()-t0:.0f}s)", flush=True)
    if not used:
        return
    t1 = time.time()
    Gk = [[cols[j].get(e, 0) for j in range(len(used))] for e in KEEP]
    ker = int_kernel(Gk) if Gk else [[1 if i == j else 0 for j in range(len(used))]
                                     for i in range(len(used))]
    print(f"  integer kernel of the KEEP block: dim {len(ker)} ({time.time()-t1:.0f}s)", flush=True)
    if not ker:
        return
    # Stage A: inside the kernel, make the failing values 0 mod p
    Gf = [[sum(cols[j].get(e, 0) * ker[t][j] for j in range(len(used))) % P
           for t in range(len(ker))] for e in FAIL]
    bf = [eqs(e) for e in FAIL]
    y = gf_solve(Gf, [(-b) % P for b in bf], P)
    print(f"  Stage A (0 mod p inside the kernel): {'SOLVABLE' if y is not None else 'not solvable'}")
    if y is None:
        return
    kk = [sum(ker[t][j] * y[t] for t in range(len(ker))) for j in range(len(used))]
    snap = [v[u] for u in used]
    for j, u in enumerate(used):
        v[u] += kk[j]
    AV2 = [atomval(a, v) for a in range(L.NA)]
    def eqs2(e): return sum(c * AV2[a] for a, c in L.eq_atoms[e][2].items())
    F2 = [e for e in range(L.NEQ) if eqs2(e) != 0]
    pdiv = sum(1 for e in F2 if eqs2(e) % P == 0)
    print(f"  after Stage A: failing={len(F2)} score={L.NEQ-len(F2)} ; "
          f"{pdiv} of {len(F2)} failing values are multiples of p", flush=True)
    if len(F2) > len(FAIL) * 3:
        print("  (kernel step was not exact -- the map is polynomial; reverting)")
        for j, u in enumerate(used):
            v[u] = snap[j]
        return
    # Stage B: exact integer solve now that the rhs is p-divisible
    E2, idx2, cols2, used2 = build_linear(v, AV2, F2)
    M = [[cols2[j].get(e, 0) for j in range(len(cols2))] for e in E2]
    rhs = [-eqs2(e) for e in E2]
    print(f"  Stage B system {len(E2)} x {len(cols2)}", flush=True)
    t2 = time.time()
    x = solve_int(M, rhs)
    print(f"  Stage B exact integer solve: {'FOUND' if x else 'none'} ({time.time()-t2:.0f}s)")
    if x:
        for j, u in enumerate(used2):
            v[u] += x[j]
        AV3 = [atomval(a, v) for a in range(L.NA)]
        F3 = [e for e in range(L.NEQ)
              if sum(c * AV3[a] for a, c in L.eq_atoms[e][2].items()) != 0]
        print(f"  FINAL failing={len(F3)} score={L.NEQ-len(F3)}")
        if len(F3) < len(FAIL):
            json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                      open(os.path.join(HERE, 'data', f'ip17_{name.replace(" ", "_")}.json'), 'w'))
            print("  saved")


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    for rel, nm in [(os.path.join(HERE, 'data', 'finish3_named.json'), 's11best39018'),
                    (os.path.join(LAB, 'best', 'new_instance_partial_39026.json'), 'checkpoint39026'),
                    (os.path.join(HERE, 'data', 'closehit2.json'), 'closehit2')]:
        try:
            run(rel, nm)
        except Exception as e:
            print(nm, 'error', repr(e))
