#!/usr/bin/env python3
"""check_linkA.py -- machine-checkable CERTIFICATE for LINK A of the chain

        QUBO ground state (E=0)   <=>   all arithmetic gadget constraints hold.

This file is READ-ONLY with respect to every existing module: it imports
qubo.py, ladder.py, squeeze/mmqb.py, squeeze/mm.py, squeeze/verify.py and never
edits them.  It monkeypatches ONE method (QB.add_square) at runtime, and only to
*record* the squares the compiler already builds -- the numerical behaviour of
the compiler is untouched (the patch calls the original and returns its result).

What is certified, per built instance Q (after finalize):

  (D) DECOMPOSITION IDENTITY (algebraic, exact, all sizes)
        Q(x)  ==  sum_k  L_k(x)^2               (square gadgets, coef >= 0)
                + W * sum_g R_g(x)              (AND/Rosenberg gadgets, W >= 1)
                + sum_o  T_o(x)                 (thermometer/order gadgets, unary)
      verified by re-expanding every recorded gadget and checking the resulting
      monomial->coefficient dict is IDENTICAL to Q.Q.  Because every summand is
      >= 0 for binary x and W >= 1, this identity is a *proof* that
        E(x)=0  <=>  every L_k(x)=0 and every R_g(x)=0 and every T_o(x)=0
      i.e. E=0 exactly on the gadget-constraint set, for ALL x (any size).

  (S) STRUCTURAL AUDIT
        - every recorded penalty is one of the three certified gadget types
          with the right coefficients (Rosenberg coefficients checked against
          the canonical a*b-2az-2bz+3z),
        - AND outputs are never AND inputs (no nested ANDs) -- the assumption
          the finalize() W_and bound silently relies on,
        - W_and audit: recompute the local-load bound and confirm
          W > max_v load(v) with an explicit margin (the rigidity proof).

  (E) EXHAUSTIVE ENUMERATION (independent oracle, small instances)
        - full 2^n brute force enumerates {x : E(x)=0} directly from Q.Q,
        - a constraint-DFS independently enumerates {x : gadgets hold},
        - assert the two sets are EQUAL.
      This confirms (D)'s conclusion with no algebra, by literal enumeration.

  (F) FAITHFULNESS (Link A -> problem), cross-checked against squeeze/verify.py
        - project the gadget-constraint set onto (a,b,c) and confirm it equals
          { (a,b,c) : a*b == c (mod p) }, and that it agrees with verify.L0X.

Run:  python3 check_linkA.py
"""
import itertools
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ANNEAL = os.path.dirname(os.path.dirname(HERE))          # .../anneal
SQUEEZE = os.path.join(ANNEAL, 'squeeze')
for d in (ANNEAL, SQUEEZE):
    if d not in sys.path:
        sys.path.insert(0, d)

import qubo                                              # noqa: E402
from qubo import QB                                      # noqa: E402

sys.setrecursionlimit(1_000_000)


# --------------------------------------------------------------------------
# Runtime recording of the squares the compiler builds (no file is modified).
# Base QB has no squares list; MMQB already keeps self.squares.  We patch only
# to record, and only for *exact* QB instances, so MMQB is never double-counted.
# --------------------------------------------------------------------------
_orig_add_square = QB.add_square


def _recording_add_square(self, lin, const, tgt=None):
    if tgt is None and type(self) is QB:
        if not hasattr(self, '_cert_squares'):
            self._cert_squares = []
        self._cert_squares.append((dict(lin), const))
    return _orig_add_square(self, lin, const, tgt)


QB.add_square = _recording_add_square


def get_squares(Q):
    """the list of (lin, const) square gadgets built into Q."""
    s = getattr(Q, 'squares', None)
    if s:
        return s
    return getattr(Q, '_cert_squares', [])


# --------------------------------------------------------------------------
# Gadget expansions (as monomial -> integer coefficient dicts over binary vars)
# --------------------------------------------------------------------------
def expand_square(lin, const):
    """(sum_v lin[v] v + const)^2  with v^2 = v  ->  monomial dict.
       Mirrors QB.add_square exactly."""
    out = defaultdict(int)
    out[()] += const * const
    vs = sorted(lin)
    for v in vs:
        c = lin[v]
        out[(v,)] += c * c + 2 * const * c
    for a in range(len(vs)):
        for b in range(a + 1, len(vs)):
            out[(vs[a], vs[b])] += 2 * lin[vs[a]] * lin[vs[b]]
    return out


def rosenberg(i, j, z):
    """canonical AND penalty  i*j - 2 i z - 2 j z + 3 z  ->  monomial dict."""
    out = defaultdict(int)
    key = (i, j) if i < j else (j, i)
    out[key] += 1
    out[(min(i, z), max(i, z))] += -2
    out[(min(j, z), max(j, z))] += -2
    out[(z,)] += 3
    return out


def prune(d):
    return {m: c for m, c in d.items() if c}


# --------------------------------------------------------------------------
# (D) decomposition identity
# --------------------------------------------------------------------------
def gadget_reconstruct(Q):
    """re-expand every gadget and sum, with the AND weight W the compiler chose."""
    R = defaultdict(int)
    for lin, const in get_squares(Q):
        for m, c in expand_square(lin, const).items():
            R[m] += c
    for (i, j), z in Q.andcache.items():
        for m, c in rosenberg(i, j, z).items():
            R[m] += Q.W * c
    for (u, w) in getattr(Q, 'orders', []):          # unary thermometer: w*(1-u)
        R[(w,)] += 1
        R[(min(u, w), max(u, w))] += -1
    return prune(R)


def check_decomposition(Q):
    R = gadget_reconstruct(Q)
    QQ = prune(dict(Q.Q))
    return R == QQ, R, QQ


# --------------------------------------------------------------------------
# (S) structural audit
# --------------------------------------------------------------------------
def check_and_coeffs(Q):
    """every AND penalty has exactly the canonical Rosenberg coefficients."""
    A = defaultdict(int)
    for (i, j), z in Q.andcache.items():
        for m, c in rosenberg(i, j, z).items():
            A[m] += c
    return prune(A) == prune(dict(Q.andpen))


def check_no_nested_ands(Q):
    outs = set(Q.andcache.values())
    ins = set()
    for (i, j) in Q.andcache:
        ins.add(i)
        ins.add(j)
    return ins & outs                                # empty set == good


def check_squares_integral(Q):
    for lin, const in get_squares(Q):
        if not isinstance(const, int):
            return False
        for v, c in lin.items():
            if not isinstance(c, int) or not isinstance(v, int):
                return False
    return True


def wand_audit(Q):
    """recompute the finalize() local-load bound and the rigidity margin.

    Flipping AND var z (its inputs a,b fixed, all else fixed) moves its own
    Rosenberg penalty by >= 1 (min of |3-2a-2b| over a,b in {0,1}), so the
    energy changes by >= W from the AND part; the non-AND (square) part changes
    by at most load(z) = sum of |coef| of pen-monomials containing z.  Because
    no AND output is an AND input (checked separately), z occurs in NO other
    Rosenberg penalty, so load(z) is the *complete* competing sensitivity.
    Hence W > max_z load(z) makes z=ab forced in every ground state."""
    andvars = set(Q.andcache.values())
    load = defaultdict(int)
    for m, c in Q.pen.items():
        for v in m:
            if v in andvars:
                load[v] += abs(c)
    maxload = max(load.values()) if load else 0
    return Q.W, maxload, Q.W - maxload               # W, maxload, margin


# --------------------------------------------------------------------------
# (E) enumeration oracles
# --------------------------------------------------------------------------
def brute_zero_states(Q):
    """{x : E(x)=0} by literal 2^n enumeration of the finalized Hamiltonian."""
    n = Q.n
    off = Q.Q.get((), 0)
    lin = [(m[0], c) for m, c in Q.Q.items() if len(m) == 1]
    qua = [(m[0], m[1], c) for m, c in Q.Q.items() if len(m) == 2]
    out = set()
    for bits in range(1 << n):
        e = off
        for v, c in lin:
            if (bits >> v) & 1:
                e += c
        for u, v, c in qua:
            if (bits >> u) & 1 and (bits >> v) & 1:
                e += c
        if e == 0:
            out.add(tuple((bits >> v) & 1 for v in range(n)))
    return out


def constraint_states(Q):
    """{x : every square = 0, every AND correct, every order ok} by DFS with
    interval propagation.  Enumerates the solution set in time ~ its size, and
    is completely independent of Q.Q (uses only the recorded gadgets)."""
    eqs = get_squares(Q)
    gates = [(z, i, j) for (i, j), z in Q.andcache.items()]
    orders = getattr(Q, 'orders', [])
    n = Q.n
    occ = [[] for _ in range(n)]
    for e, (lin, k) in enumerate(eqs):
        for v in lin:
            occ[v].append(e)
    val = [None] * n
    out = []

    def feasible(elist):
        for e in elist:
            lin, k = eqs[e]
            lo = hi = k
            for v, c in lin.items():
                if val[v] is None:
                    lo += min(0, c)
                    hi += max(0, c)
                else:
                    lo += c * val[v]
                    hi += c * val[v]
            if lo > 0 or hi < 0:
                return False
        return True

    def rec(v):
        if v == n:
            out.append(tuple(val))
            return
        for b in (0, 1):
            val[v] = b
            ok = feasible(occ[v])
            if ok:
                for (z, i, j) in gates:
                    a, bb, cc = val[i], val[j], val[z]
                    if a is not None and bb is not None and cc is not None and cc != a * bb:
                        ok = False
                        break
                    if cc == 1 and ((a == 0) or (bb == 0)):
                        ok = False
                        break
            if ok:
                for (u, w) in orders:
                    if val[u] == 0 and val[w] == 1:
                        ok = False
                        break
            if ok:
                rec(v + 1)
        val[v] = None

    rec(0)
    return set(out)


# --------------------------------------------------------------------------
# per-instance driver
# --------------------------------------------------------------------------
class Fail(Exception):
    pass


def certify_instance(Q, label, brute_cap=21, enumerate_constraints=True):
    """Run (D),(S) always; (E) full-brute if n<=brute_cap; return a report line.

    (D) is a mathematical identity that proves {E=0}=={constraints} at ANY size,
    so the constraint enumeration is only extra confidence.  For large instances
    the constraint set carries range-slack ancilla multiplicity and can be huge;
    pass enumerate_constraints=False there and rely on (D)."""
    # (S) structural
    if not check_squares_integral(Q):
        raise Fail(f"{label}: a square gadget has non-integer coefficients")
    if not check_and_coeffs(Q):
        raise Fail(f"{label}: an AND penalty is not the canonical Rosenberg form")
    nest = check_no_nested_ands(Q)
    if nest:
        raise Fail(f"{label}: AND output(s) used as AND input(s): {nest}")
    # (D) decomposition identity
    ok, R, QQ = check_decomposition(Q)
    if not ok:
        diff = {m: (R.get(m, 0), QQ.get(m, 0)) for m in set(R) | set(QQ)
                if R.get(m, 0) != QQ.get(m, 0)}
        raise Fail(f"{label}: decomposition != Q.Q; first diffs "
                   f"{list(diff.items())[:6]}")
    # (S) W_and audit
    W, maxload, margin = wand_audit(Q)
    if margin < 1:
        raise Fail(f"{label}: W_and={W} <= max local load {maxload}; NOT RIGID")
    csts = None
    if enumerate_constraints:
        # constraint enumeration + forward check every state has E=0
        csts = constraint_states(Q)
        for x in csts:
            if Q.energy(list(x)) != 0:
                raise Fail(f"{label}: a gadget-constraint state has E != 0")
        # (E) independent full brute for small n
        if Q.n <= brute_cap:
            zs = brute_zero_states(Q)
            if zs != csts:
                raise Fail(f"{label}: {{E=0}} != {{constraints}} "
                           f"(|E=0|={len(zs)} |constraints|={len(csts)})")
            tag = f"E ok (2^{Q.n} brute, |set|={len(csts)})"
        else:
            tag = (f"D-certified; forward-checked |constraints|={len(csts)}")
    else:
        tag = "D+S-certified (identity proves {E=0}=={constraints})"
    return csts, (f"  {label:52s} n={Q.n:4d} squares={len(get_squares(Q)):4d} "
                  f"ands={len(Q.andcache):3d} W={W:<4d} margin={margin:<4d} {tag}")


# --------------------------------------------------------------------------
# gadget UNIT certifications (each of the three types, in isolation)
# --------------------------------------------------------------------------
def unit_square():
    """(2a + 3b + c - 1)^2 >= 0, == 0 iff 2a+3b+c=1.  Full 8-state brute."""
    Q = QB(mode='binary')
    a, b, c = Q.new('a', 'x'), Q.new('b', 'x'), Q.new('c', 'x')
    Q.add_square({a: 2, b: 3, c: 1}, -1)
    Q.finalize()
    bad = 0
    for va, vb, vc in itertools.product((0, 1), repeat=3):
        e = Q.energy([va, vb, vc])
        want0 = (2 * va + 3 * vb + vc - 1) == 0
        if (e == 0) != want0 or e < 0:
            bad += 1
    return bad == 0


def unit_and():
    """z = AND(a,b): W*(ab-2az-2bz+3z) >= 0, ==0 iff z=ab.  Full brute."""
    Q = QB(mode='binary')
    a, b = Q.new('a', 'x'), Q.new('b', 'x')
    z = Q.AND(a, b)
    Q.finalize()
    bad = 0
    for va, vb, vz in itertools.product((0, 1), repeat=3):
        x = [0, 0, 0]
        x[a], x[b], x[z] = va, vb, vz
        e = Q.energy(x)
        if (e == 0) != (vz == va * vb) or e < 0:
            bad += 1
    return bad == 0


def unit_onehot(D):
    """Exactly-one over D selectors via the build_win sequential-counter chain.
       Full brute over all 2^(2D-1) variables; project E=0 onto the selectors
       and confirm it equals the set of one-hot vectors."""
    Q = QB(mode='binary')
    u = [Q.new(f"u{t}", 'input') for t in range(D)]
    prev = None
    for t, v in enumerate(u[:-1]):
        pv = Q.new(f"p{t}", 'onehot')
        lin = {pv: -1, v: 1}
        if prev is not None:
            lin[prev] = 1
        Q.add_square(lin, 0)
        prev = pv
    Q.add_square({prev: 1, u[-1]: 1}, -1)
    Q.finalize()
    n = Q.n
    sel_zero = set()
    for bits in range(1 << n):
        x = [(bits >> v) & 1 for v in range(n)]
        if Q.energy(x) == 0:
            sel_zero.add(tuple(x[v] for v in u))
    onehot = {tuple(1 if i == k else 0 for i in range(D)) for k in range(D)}
    return sel_zero == onehot


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    import verify                                       # squeeze/verify.py

    print("=" * 92)
    print("LINK A CERTIFICATE  --  QUBO E=0  <=>  gadget constraints")
    print("=" * 92)

    print("\n[1] GADGET UNIT CERTIFICATIONS (each type in isolation, full brute)")
    us, ua = unit_square(), unit_and()
    print(f"  (a) add_square  (2a+3b+c-1)^2 >= 0, =0 iff form=0 ....... "
          f"{'OK' if us else '*** FAIL ***'}")
    print(f"  (b) AND penalty  z=ab, W*(ab-2az-2bz+3z) >= 0 .......... "
          f"{'OK' if ua else '*** FAIL ***'}")
    oh_ok = True
    for D in (2, 3, 4, 5, 6, 8):
        r = unit_onehot(D)
        oh_ok = oh_ok and r
        print(f"  (c) one-hot/seq-counter  D={D}: {{E=0}}|_sel == one-hot . "
              f"{'OK' if r else '*** FAIL ***'}")
    fails = 0 if (us and ua and oh_ok) else 1

    print("\n[2] MODMUL INSTANCES  --  decomposition identity, W_and audit,")
    print("    and {E=0} == {constraints} by independent 2^n brute where feasible.")
    modmul_specs = []
    for p in (2, 3):                                    # small enough for full brute
        for mode in ('binary', 'wallace'):
            for red in ('naf', 'quotient'):
                modmul_specs.append((p, 'schoolbook', red, mode))
    for p in (5, 7, 13):                                # decomposition-certified
        for mode in ('binary', 'wallace'):
            for mult in ('schoolbook', 'karatsuba', 'toom3'):
                modmul_specs.append((p, mult, 'naf', mode))
    for p in (29, 61, 127, 251):                        # larger: identity + forward
        modmul_specs.append((p, 'schoolbook', 'naf', 'wallace'))

    for (p, mult, red, mode) in modmul_specs:
        Q, A, B, C, base = verify.make(p, mult=mult, leaf=3, red=red, mode=mode)
        label = f"p={p} {mult} {red} {mode}"
        # enumerate the constraint set only when it is cheap (small p); for
        # larger p the decomposition identity (D) already proves set equality.
        enum = p <= 13
        try:
            csts, line = certify_instance(Q, label, enumerate_constraints=enum)
            print(line)
        except Fail as ex:
            print(f"  *** FAIL *** {ex}")
            fails += 1

    print("\n[2b] LARGER p, independent input-exhaustive cross-check via verify.L1")
    print("     (every (a,b); correct c plus wrong c; no spurious E=0 admitted).")
    for p in (29, 61):
        for mode in ('binary', 'wallace'):
            _, ok, bad = verify.L1(p, mult='schoolbook', leaf=3, red='naf',
                                   mode=mode)
            if bad:
                fails += 1
            print(f"  p={p:3d} schoolbook naf {mode:7s}  checked={ok + bad:6d} "
                  f"bad={bad}  {'OK' if bad == 0 else '*** FAIL ***'}")
    for p in (127, 251):
        _, ok, bad = verify.L1(p, sample=8, mult='schoolbook', leaf=3, red='naf',
                               mode='wallace')
        if bad:
            fails += 1
        print(f"  p={p:3d} schoolbook naf wallace  checked={ok + bad:6d} "
              f"bad={bad}  {'OK' if bad == 0 else '*** FAIL ***'}")

    print("\n[3] FAITHFULNESS  --  {constraints}|_(a,b,c) == {a*b==c (mod p)},")
    print("    cross-checked against verify.L0X.")
    for p in (3, 5, 7, 13):
        for mult in ('schoolbook', 'karatsuba', 'toom3'):
            for mode in ('binary', 'wallace'):
                Q, A, B, C, base = verify.make(p, mult=mult, leaf=3, red='naf',
                                               mode=mode)
                csts = constraint_states(Q)
                s = p.bit_length()
                proj = set()
                for x in csts:
                    a = sum(x[v] << t for t, v in enumerate(A.bits))
                    b = sum(x[v] << t for t, v in enumerate(B.bits))
                    c = sum(x[v] << t for t, v in enumerate(C.bits))
                    proj.add((a, b, c))
                truth = {(a, b, c)
                         for a in range(1 << s) for b in range(1 << s)
                         for c in range(1 << s) if (a * b - c) % p == 0}
                # cross-check against verify's own enumeration
                _, vzeros, vtruth = verify.L0X(p, mult=mult, leaf=3, red='naf',
                                               mode=mode)
                ok = (proj == truth) and (proj == vzeros) and (truth == vtruth)
                if not ok:
                    fails += 1
                print(f"  p={p:2d} {mult:10s} {mode:7s}  |constraints|-proj={len(proj):4d}"
                      f"  |truth|={len(truth):4d}  verify-agrees={proj == vzeros}  "
                      f"{'FAITHFUL' if ok else '*** MISMATCH ***'}")

    print("\n[4] LADDER ONE-HOT WINDOW INSTANCE (base QB, sequential-counter gadget)")
    try:
        from ecsmall import curve, find
        from ladder import build_win
        p, B, m, w = 97, 3, 4, 2
        add, mul = curve(p, B)
        G, order = find(p, B)
        M = (m + w - 1) // w
        table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)]
                 for j in range(M)]
        off = sum(1 << (w * j) for j in range(M))
        # pick a k whose window chain is non-degenerate
        def chain_ok(dg):
            S = table[0][dg[0]]
            for j in range(1, M):
                Qp = table[j][dg[j]]
                if S is None or Qp is None or S[0] == Qp[0]:
                    return False
                S = add(S, Qp)
            return S is not None
        dgs = lambda kk: [(kk >> (w * j)) % (1 << w) for j in range(M)]
        k = next(kk for kk in range(1 << m) if chain_ok(dgs(kk)))
        Tp = add(mul(k % order, G), mul(off % order, G))
        L, U = build_win(p, B, table, Tp, w, mode='wallace')
        Ql = L.qb
        # structural + decomposition + W audit + forward (n too big for full brute)
        if not check_and_coeffs(Ql):
            raise Fail("ladder: AND penalty not canonical")
        if check_no_nested_ands(Ql):
            raise Fail("ladder: nested ANDs present")
        ok, R, QQ = check_decomposition(Ql)
        if not ok:
            raise Fail("ladder: decomposition != Q.Q")
        W, maxload, margin = wand_audit(Ql)
        if margin < 1:
            raise Fail(f"ladder: W_and margin {margin} < 1")
        # witness for the true k must be E=0 (forward direction sanity)
        dg = dgs(k)
        wv0 = {f"_u{j}": dg[j] for j in range(M)}
        inp = {}
        for j in range(M):
            for t in range(1 << w):
                inp[U[j][t]] = 1 if t == dg[j] else 0
        x, _ = Ql.witness(inp, wv0)
        assert Ql.energy(x) == 0, "true-k witness not E=0"
        # count one-hot squares that are present (sequential counter chain)
        onehot_sq = sum(1 for lin, c in get_squares(Ql)
                        if c in (0, -1) and len(lin) <= 3)
        print(f"  build_win p={p} m={m} w={w}: n={Ql.n} squares={len(get_squares(Ql))} "
              f"ands={len(Ql.andcache)} W={W} margin={margin}")
        print(f"    decomposition identity Q == sum(squares)+W*sum(AND): OK")
        print(f"    W_and rigidity  W > max local load ({maxload}): OK")
        print(f"    true-scalar witness energy: {Ql.energy(x)}  (E=0): OK")
    except Fail as ex:
        print(f"  *** FAIL *** {ex}")
        fails += 1
    except Exception as ex:
        print(f"  *** ERROR *** {type(ex).__name__}: {ex}")
        fails += 1

    print("\n" + "=" * 92)
    print(f"TOTAL FAILURES: {fails}")
    print("=" * 92)
    return fails


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
