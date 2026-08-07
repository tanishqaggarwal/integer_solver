"""LATTICE-COLUMN PRICER — knobs chosen by which lattice column they move.

Why this instrument exists
--------------------------
Every pricer I have built so far chooses knobs by HANDLE INCIDENCE: a subset W names p-handles,
`ieng.site` frees them plus their definer-level collateral, and those five-or-so variables are
the only columns the solver ever sees.  The 2^12 / 2^16 / 2^18 enumerations are complete under
four such instruments and nothing exceeds 39,026.

Agent W measured why that space cannot contain the answer: inside K the system is CONSISTENT
over Q -- `rank[A|b] = 28` with the rhs not a pivot -- and blocked only over Z.  **The binding
object is a lattice index, not a rank.**  A lattice index is moved by adding columns that change
the elementary divisors, and that is a property of a column, not of a handle.  Handle incidence
cannot select for it, so the knob set was the wrong object.

What this module does differently
---------------------------------
1. **Frame D**: it linearises at the DELIVERABLE itself (39,026, 7 failing), not at the
   uncorrupted baseline.  Rows are the 7 failing equations; the rhs is the deliverable's own
   residual.  That is the system in which "buy a failing equation" is the question.
2. **Candidate pool by reachability, not incidence**: every free input that can reach an atom of
   a target equation, found by ONE backward walk over the definer graph.
3. **Selection by lattice contribution**: a column is scored by what it does to the Smith normal
   form of the target system -- does it raise the rank, or divide down an elementary divisor and
   so shrink the index?  Columns that change nothing about the lattice are discarded no matter
   how incident they look.
4. **Price by re-propagation, as always.**  Collateral damage is measured, never assumed.

Exactness discipline is inherited unchanged: the columns are exact integer derivatives obtained
by re-propagating the real engine and verified affine by second differences; anything that fails
the affinity test is dropped rather than approximated.
"""
import sys, os, json, time, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H                                            # noqa: E402
import ieng, fscore, sparse                                    # noqa: E402
import flint                                                   # noqa: E402

NV = ieng.NV
D4 = [642, 28730, 29854, 31864]
DELIV_FAILS = [12231, 12270, 12350, 14584, 18673, 22044, 29125]

# ---------------------------------------------------------------- the two frames
FREED, PIN = ieng.site(D4)
VD = ieng.VD

_sc, BAD_D, V_D = ieng.score_from_unc({u: VD[u] for u in FREED}, PIN)
SCORE_D = _sc
FAILS_D = sorted(fscore.fails(BAD_D))


def eqmaps(fails):
    """atom -> coefficient, plus the constant part, for each equation in `fails`."""
    out = {}
    for e in fails:
        cm = collections.defaultdict(int); c0 = 0
        for c, a in H.eqt[e][2]:
            if a < 0:
                c0 += c
            else:
                cm[a] += c
        out[e] = (dict(cm), c0)
    return out


CM_D = eqmaps(FAILS_D)


def eqval(bad, e, cmd=None):
    """exact integer value of equation e given the nonzero-atom map `bad`."""
    cm, c0 = (cmd or CM_D)[e]
    return c0 + sum(c * bad[a] for a, c in cm.items() if a in bad)


# ---------------------------------------------------------------- candidate pool
def reach_pool(fails):
    """Every FREE INPUT that can reach an atom of one of `fails`, by one backward walk.

    Edge: w's definer atom mentions u  =>  u influences w.  So from the variables of the
    target atoms, walk back through definer atoms; every source reached is a free input.
    This is REACHABILITY (can this input move the row at all), which is a necessary
    condition and cheap; the lattice criterion below is what actually selects.
    """
    atoms = set()
    for e in fails:
        for _c, a in H.eqt[e][2]:
            if a >= 0:
                atoms.add(a)
    seen = set(); st = []
    for a in atoms:
        for u in H.avars[a]:
            if u not in seen:
                seen.add(u); st.append(u)
    free = set()
    while st:
        w = st.pop()
        d = H.definer[w]
        if d is None or w in PIN:
            free.add(w); continue
        i, _k = d
        for u in H.avars[i]:
            if u != w and u not in seen:
                seen.add(u); st.append(u)
    return sorted(free), atoms


# ---------------------------------------------------------------- exact columns
def column(u, base_v, base_bad, fails, cmd, pin, probes=(1, 2, 7)):
    """Exact integer derivative of each target equation w.r.t. free input u.

    Affinity is TESTED, not assumed: the atom deltas must be exactly linear in the probe at
    +1, +2 and +7.  Returns None if any atom fails, so a nonlinear input is dropped rather
    than linearised badly.
    """
    o = base_v[u]
    try:
        bs = [ieng.resid(base_v, base_bad, {u: o + p}, pin)[0] for p in probes]
    except Exception:
        return None
    d1 = {}
    keys = set(base_bad)
    for b in bs:
        keys |= set(b)
    for a in keys:
        z = base_bad.get(a, 0)
        v1 = bs[0].get(a, 0) - z
        for k, p in enumerate(probes[1:], start=1):
            if bs[k].get(a, 0) - z != p * v1:
                return None
        if v1:
            d1[a] = v1
    col = {}
    for e in fails:
        cm, _c0 = cmd[e]
        s = 0
        for a, dv in d1.items():
            c = cm.get(a)
            if c:
                s += c * dv
        if s:
            col[e] = s
    return col


# ---------------------------------------------------------------- lattice algebra
def _mat(cols, knobs, fails):
    return flint.fmpz_mat([[cols[u].get(e, 0) for u in knobs] for e in fails])


def snf_report(cols, knobs, rhs, fails):
    """Everything the lattice says about A d = b, computed exactly over Z.

    Returns rank(A), rank([A|b]) -- so the Q-consistency question W answered is answered here
    in MY parse -- the elementary divisors of A, the index (their product), and whether the
    system is solvable over Z.
    """
    A = _mat(cols, knobs, fails)
    Ab = flint.fmpz_mat([[cols[u].get(e, 0) for u in knobs] + [rhs[e]] for e in fails])
    rA, rAb = A.rank(), Ab.rank()
    S = A.snf()
    divs = [S[i, i] for i in range(min(S.nrows(), S.ncols())) if S[i, i] != 0]
    idx = 1
    for d in divs:
        idx *= int(d)
    solvable = (rA == rAb) and _int_solvable(cols, knobs, rhs, fails)
    return {'nrows': len(fails), 'ncols': len(knobs), 'rank_A': rA, 'rank_Ab': rAb,
            'rhs_is_pivot': rAb > rA, 'divisors': [int(d) for d in divs],
            'index': idx, 'Z_solvable': solvable}


def _int_solvable(cols, knobs, rhs, fails):
    """A d = b over Z, decided by the standard SNF criterion but computed via HNF of [A|b].

    d(k) invariants: A d = b is solvable over Z iff for every k, gcd of k x k minors of A
    equals that of [A|b] (and the ranks agree).  Comparing the products of the nonzero
    elementary divisors of A and [A|b] IS that test, since SNF divisors are ratios of the d(k).
    """
    A = _mat(cols, knobs, fails)
    Ab = flint.fmpz_mat([[cols[u].get(e, 0) for u in knobs] + [rhs[e]] for e in fails])
    if A.rank() != Ab.rank():
        return False
    da = [int(A.snf()[i, i]) for i in range(min(A.nrows(), A.ncols()))]
    db = [int(Ab.snf()[i, i]) for i in range(min(Ab.nrows(), Ab.ncols()))]
    da = [x for x in da if x] + [0] * 99
    db = [x for x in db if x] + [0] * 99
    return all(da[i] == db[i] for i in range(A.rank()))


def index_of(cols, knobs, fails):
    if not knobs:
        return None
    A = _mat(cols, knobs, fails)
    S = A.snf()
    idx = 1
    for i in range(min(S.nrows(), S.ncols())):
        d = int(S[i, i])
        if d:
            idx *= d
    return A.rank(), idx


def select_by_lattice(cols, pool, rhs, fails, cap=60, verbose=True):
    """Greedy knob selection by LATTICE CONTRIBUTION, not by incidence.

    A candidate is taken when it raises the rank of A, or (rank unchanged) when it strictly
    shrinks the index -- the product of the elementary divisors, which is exactly the
    obstruction W identified.  A column that leaves both alone is discarded however many
    target rows it touches.  Stops as soon as the system is solvable over Z.
    """
    knobs = []
    cur = (0, 0)
    for step in range(cap):
        best = None
        for u in pool:
            if u in knobs or not cols[u]:
                continue
            r, i = index_of(cols, knobs + [u], fails)
            if r > cur[0] or (r == cur[0] and cur[1] and i < cur[1]):
                key = (r, -i)
                if best is None or key > best[0]:
                    best = (key, u, (r, i))
        if best is None:
            break
        knobs.append(best[1]); cur = best[2]
        if verbose:
            print(f'  +x_{best[1]:<6d} rank {cur[0]}  index {cur[1]}', flush=True)
        if _int_solvable(cols, knobs, rhs, fails):
            if verbose:
                print('  -> integrally solvable', flush=True)
            break
    return knobs


# ---------------------------------------------------------------- price
def price(knobs, cols, rhs, fails, base_v, base_bad, pin, want=False):
    """Solve A d = b over Z on `knobs`, APPLY it, re-propagate, and measure the real score."""
    rows = [{u: cols[u][e] for u in knobs if e in cols[u]} for e in fails]
    rr = [rhs[e] for e in fails]
    keep, sols = [], []
    for i in range(len(rows)):
        if not rows[i]:
            continue
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rr[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial; sols.append(s)
    out = {'nrows_solved': len(keep), 'best': None}
    best = (fscore.score(base_bad), None)
    for s in sols:
        ch = {u: base_v[u] + d for u, d in s.items() if d}
        if not ch:
            continue
        try:
            bad, v = ieng.resid(base_v, base_bad, ch, pin)
            sc = fscore.score(bad)
        except Exception:
            continue
        if sc > best[0]:
            best = (sc, ch)
    out['score'] = best[0]
    if want:
        out['changes'] = best[1]
    return out
