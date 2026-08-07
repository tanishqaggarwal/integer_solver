"""IP #3 -- MINIMUM-WEIGHT COSET.

At a state the reachable perturbations form an affine integer lattice.  Writing b for the
vector of equation inner-sums and G for the matrix of handle effects on them, the objective is
exactly

        minimise   || b + G k ||_0     over integer k

i.e. a minimum-weight coset problem.  Constrain the currently-satisfied equations to stay zero
(that is the integer kernel of G restricted to them), then minimise the number of nonzeros on
the failing ones inside that kernel.  Exact, by HNF + subset search.
"""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
sys.set_int_max_str_digits(400000)


def int_kernel(M):
    """integer kernel basis of M (rows = equations, cols = variables) via column HNF."""
    m = len(M)
    n = len(M[0]) if m else 0
    A = [row[:] for row in M]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    r = 0
    for row in range(m):
        while True:
            nz = [j for j in range(r, n) if A[row][j] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda j: abs(A[row][j]))
            j0 = nz[0]
            moved = False
            for j in nz[1:]:
                q = A[row][j] // A[row][j0]
                if q:
                    for i in range(m):
                        A[i][j] -= q * A[i][j0]
                    for i in range(n):
                        U[i][j] -= q * U[i][j0]
                    moved = True
            if not moved:
                break
        nz = [j for j in range(r, n) if A[row][j] != 0]
        if not nz:
            continue
        j0 = nz[0]
        if j0 != r:
            for i in range(m):
                A[i][r], A[i][j0] = A[i][j0], A[i][r]
            for i in range(n):
                U[i][r], U[i][j0] = U[i][j0], U[i][r]
        r += 1
    return [[U[i][j] for i in range(n)] for j in range(r, n)]   # kernel basis vectors


def load(name):
    v = [0] * L.NVARS
    for k, x in json.load(open(os.path.join(HERE, 'data', name))).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    fw.forward(v)
    return v


def run(name, maxfail=10):
    v = load(name)
    av = L.all_atom_values(v)
    FAIL = sorted(L.failing_eqs(av))
    print(f"=== {name}: failing={len(FAIL)}  score={L.NEQ-len(FAIL)}")

    # candidate handles: every free input touching any atom of any equation we could disturb
    region_atoms = set()
    for e in FAIL:
        region_atoms |= set(L.eq_atoms[e][2])
    cands = set()
    for a in region_atoms:
        for u in L.avars[a]:
            if L.definer.get(u) is None:
                cands.add(u)
    cands = sorted(cands)

    # equations any candidate can touch
    E = set(FAIL)
    for t in cands:
        for a in L.var_atoms[t]:
            E |= set(L.atom2eq.get(a, {}))
    E = sorted(E)

    def inner(vv):
        a2 = L.all_atom_values(vv)
        return [sum(c * a2[a] for a, c in L.eq_atoms[e][2].items()) for e in E]

    b = inner(v)
    cols, used = [], []
    for t in cands:
        old = v[t]
        v[t] = old + 1
        fw.forward(v)
        e1 = inner(v)
        v[t] = old + 2
        fw.forward(v)
        e2 = inner(v)
        v[t] = old
        fw.forward(v)
        d1 = [e1[i] - b[i] for i in range(len(E))]
        d2 = [e2[i] - b[i] for i in range(len(E))]
        if any(d1) and all(d2[i] == 2 * d1[i] for i in range(len(E))):
            cols.append(d1)
            used.append(t)
    print(f"    region: {len(E)} equations, {len(used)} exact-linear handles")
    if not used:
        return len(FAIL)
    idx = {e: i for i, e in enumerate(E)}
    KEEP = [e for e in E if e not in set(FAIL)]      # must stay zero
    Gk = [[cols[j][idx[e]] for j in range(len(used))] for e in KEEP]
    ker = int_kernel(Gk) if Gk else [[1 if i == j else 0 for j in range(len(used))]
                                     for i in range(len(used))]
    print(f"    kernel dimension (moves preserving all {len(KEEP)} satisfied equations): {len(ker)}")
    if not ker:
        return len(FAIL)
    # project the failing rows onto the kernel
    Gf = [[sum(cols[j][idx[e]] * ker[t][j] for j in range(len(used))) for t in range(len(ker))]
          for e in FAIL]
    bf = [b[idx[e]] for e in FAIL]
    for allow in range(0, min(maxfail, len(FAIL)) + 1):
        for combo in itertools.combinations(range(len(FAIL)), allow):
            keep = [i for i in range(len(FAIL)) if i not in combo]
            x = solve_int([Gf[i] for i in keep], [-bf[i] for i in keep])
            if x is not None:
                print(f"    IP OPTIMUM: {allow} failing equations -> score {L.NEQ-allow}"
                      f"   (leave {[FAIL[i] for i in combo]})")
                kk = [sum(ker[t][j] * x[t] for t in range(len(ker))) for j in range(len(used))]
                for j, t in enumerate(used):
                    v[t] += kk[j]
                fw.forward(v)
                f2 = L.failing_eqs(L.all_atom_values(v))
                print(f"    applied -> failing={len(f2)} score={L.NEQ-len(f2)}")
                if len(f2) < len(FAIL):
                    json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                              open(os.path.join(HERE, 'data', f'ip3_{name}'), 'w'))
                return len(f2)
    print(f"    no improvement with <= {maxfail} allowed failures")
    return len(FAIL)


if __name__ == '__main__':
    for nm in ['finish3.json', 'closehit2.json']:
        t0 = time.time()
        run(nm)
        print(f"    ({time.time()-t0:.0f}s)", flush=True)
