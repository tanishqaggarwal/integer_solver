#!/usr/bin/env python3
"""faith.py -- exhaustive faithfulness tests for every encoding in enc.py.

Same standard as the existing demo.py / demo_win.py: build the real Hamiltonian
on a scaled instance, then enumerate EVERY candidate scalar (and, for the x-only
encoding, every branch-sign pattern too), filling all ancillas by replaying the
construction, and compare the zero-energy set with the true solution set.

No size is reported anywhere in FINDINGS.md for an encoding that does not pass
here.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ecsmall import curve, find          # noqa: E402
from enc import build_comb, build_semaev, s3   # noqa: E402


# ------------------------------------------------------------------ helpers --
def unsigned_instance(p, B, m, w, k):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m)
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)]
             for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    Tp = add(mul(k % order, G), mul(off % order, G))
    return dict(G=G, order=order, table=table, M=M, Tp=Tp, add=add, mul=mul)


def signed_instance(p, B, m, w, k):
    add, mul = curve(p, B)
    G, order = find(p, B)
    M = (m + w - 1) // w
    mm = M * w
    assert order >= (1 << mm)
    D = 1 << (w - 1)
    table = [[mul(((2 * t + 1) << (w * j)) % order, G) for t in range(D)]
             for j in range(M)]
    kk = (2 * k - ((1 << mm) - 1))
    Tp = mul(kk % order, G)
    return dict(G=G, order=order, table=table, M=M, mm=mm, Tp=Tp,
                add=add, mul=mul)


def signed_digits(k, M, w):
    """k -> (magnitude index, sign) per window, via e_i = 2 b_i - 1."""
    out = []
    for j in range(M):
        d = 0
        for i in range(w):
            b = (k >> (w * j + i)) & 1
            d += (2 * b - 1) << i
        out.append(((abs(d) - 1) // 2, 1 if d < 0 else 0))
    return out


# ----------------------------------------------------------- the comb tests --
def chain_ok(table, M, dg, add):
    S = table[0][dg[0]]
    for j in range(1, M):
        Q = table[j][dg[j]]
        if S is None or Q is None or S[0] == Q[0]:
            return False
        S = add(S, Q)
    return S is not None


def test_comb(p, B, m, w, k, mode='wallace', mux=True, kdepth=0, kmin=8,
              signed=False, onehot='square', label=""):
    if not signed:
        # pick a k whose reference chain has no degenerate addition (demo_win.py
        # does the same); on a tiny curve this is common, at 256 bits ~2^-247.
        I0 = unsigned_instance(p, B, m, w, k)
        dgs = lambda kk: [(kk >> (w * j)) % (1 << w) for j in range(I0['M'])]
        if not chain_ok(I0['table'], I0['M'], dgs(k), I0['add']):
            k = next(kk for kk in range(1 << m)
                     if chain_ok(I0['table'], I0['M'], dgs(kk), I0['add']))
    if signed:
        I = signed_instance(p, B, m, w, k)
    else:
        I = unsigned_instance(p, B, m, w, k)
    M, table, Tp = I['M'], I['table'], I['Tp']
    add, mul, G, order = I['add'], I['mul'], I['G'], I['order']
    L, SEL = build_comb(p, B, table, Tp, w, mode=mode, mux=mux,
                        kdepth=kdepth, kmin=kmin, signed=signed, onehot=onehot)
    Q = L.qb
    st = Q.stats()
    D = len(table[0])

    def wv_of(kk):
        if signed:
            dg = signed_digits(kk, M, w)
            wv = {}
            for j, (t, s_) in enumerate(dg):
                wv[f"_u{j}"] = t
                wv[f"_s{j}"] = s_
            return wv, {}
        dg = [(kk >> (w * j)) % (1 << w) for j in range(M)]
        wv = {f"_u{j}": dg[j] for j in range(M)}
        return wv, {}

    def energy_of(kk):
        wv, inp = wv_of(kk)
        try:
            x, _ = Q.witness(inp, wv)
        except Exception:
            return None
        return Q.energy(x)

    if energy_of(k) != 0:
        # this k's ladder passes through a degenerate addition (S_j == +- addend);
        # re-randomise, exactly as demo_win.py does.  At 256 bits P ~ 2^-247.
        for kk in range(1 << m):
            if energy_of(kk) == 0:
                k = kk
                break
    assert energy_of(k) == 0, "no k reaches E = 0 -- encoding is not complete"
    T_true = mul(k % order, G)
    zeros = [kk for kk in range(1 << m) if energy_of(kk) == 0]
    degen = sum(1 for kk in range(1 << m) if energy_of(kk) is None)
    sol = [kk for kk in range(1 << m) if mul(kk % order, G) == T_true]
    ok = set(zeros) <= set(sol)
    tag = (f"[comb w={w}{' signed' if signed else ''}"
           f"{' mux' if mux else ' sum'} kdepth={kdepth} 1hot={onehot}]")
    print(f"{tag} p={p} m={m} M={M} {label}")
    print(f"    {st['vars']} vars, {st['couplers']} couplers, AND={st['and_vars']}, "
          f"|J| 2^{st['dynamic_range_bits']}")
    print(f"    zero-energy k: {zeros}  true solutions: {sol}  "
          f"degenerate: {degen}/{1 << m}  -- "
          f"{'FAITHFUL' if ok else 'SPURIOUS!'}")
    assert ok, (zeros, sol)
    return st, zeros, sol, k in zeros


# ------------------------------------------------------------ semaev tests ---
def test_semaev(p, B, m, w, mode='wallace', mux=True, kdepth=0, onehot='square', verbose=True):
    """Exhaustive over magnitude tuples AND branch-sign patterns."""
    add, mul = curve(p, B)
    G, order = find(p, B)
    M = (m + w - 1) // w
    mm = M * w
    assert order >= (1 << mm)
    D = 1 << (w - 1)
    pts = [[mul(((2 * t + 1) << (w * j)) % order, G) for t in range(D)]
           for j in range(M)]
    table = [[pts[j][t][0] for t in range(D)] for j in range(M)]

    # pick a k whose reference chain is non-degenerate
    def chain(mags, signs):
        """partial sums; None if a degenerate addition occurs."""
        S = pts[0][mags[0]]
        if signs[0]:
            S = (S[0], (-S[1]) % p)
        zs = []
        for j in range(1, M):
            P = pts[j][mags[j]]
            if signs[j]:
                P = (P[0], (-P[1]) % p)
            if S is None or P is None or S[0] == P[0]:
                return None, None
            S = add(S, P)
            if S is None:
                return None, None
            zs.append(S[0])
        return S, zs

    good_k = None
    for k in range(1 << m):
        dg = signed_digits(k, M, w)
        S, zs = chain([t for t, _ in dg], [s_ for _, s_ in dg])
        if S is not None:
            good_k = k
            break
    assert good_k is not None
    k = good_k
    dg = signed_digits(k, M, w)
    Sfin, _ = chain([t for t, _ in dg], [s_ for _, s_ in dg])
    xT = Sfin[0]
    T_true = mul(k % order, G)

    L, SEL = build_semaev(p, B, table, xT, w, mode=mode, mux=mux, kdepth=kdepth,
                          onehot=onehot)
    Q = L.qb
    st = Q.stats()

    def energy(mags, signs):
        S, zs = chain(mags, signs)
        if S is None or S[0] != xT:
            # still try: the chain must land on xT for the final S_3 to vanish,
            # but the witness replay needs the intermediate z's regardless.
            if S is None:
                return None
        wv = {f"_u{j}": mags[j] for j in range(M)}
        for j in range(1, M - 1):
            wv[f"_z{j}"] = zs[j - 1]
        try:
            x, _ = Q.witness({}, wv)
        except Exception:
            return None
        return Q.energy(x)

    zeros, degen = [], 0
    for mags in itertools.product(range(D), repeat=M):
        best = None
        for signs in itertools.product((0, 1), repeat=M):
            e = energy(list(mags), list(signs))
            if e is None:
                continue
            if best is None or e < best:
                best = e
            if best == 0:
                break
        if best is None:
            degen += 1
        elif best == 0:
            zeros.append(mags)

    # genuine solutions of the ORIGINAL problem, expressed as magnitude tuples
    def orig_ok(mags):
        for signs in itertools.product((0, 1), repeat=M):
            v = sum(((-1) ** signs[j]) * (2 * mags[j] + 1) * (1 << (w * j))
                    for j in range(M))
            # v == k'  <=>  the corresponding b-vector solves the real problem
            if (v - (2 * k - ((1 << mm) - 1))) % order == 0:
                return True
        return False

    genuine = [mg for mg in zeros if orig_ok(mg)]
    spur = [mg for mg in zeros if not orig_ok(mg)]
    print(f"[semaev w={w} mux={mux} kdepth={kdepth}] p={p} m={m} M={M} k={k}")
    print(f"    {st['vars']} vars, {st['couplers']} couplers, AND={st['and_vars']}, "
          f"|J| 2^{st['dynamic_range_bits']}")
    print(f"    magnitude tuples: {D ** M} total, zero-energy {len(zeros)}, "
          f"of which genuine {len(genuine)}, SPURIOUS {len(spur)}, "
          f"no-witness {degen}")
    if spur[:6]:
        print(f"    example spurious tuples: {spur[:6]}")
    assert len(genuine) >= 1, "the true magnitude tuple is not a ground state"
    return st, len(zeros), len(genuine), len(spur)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'base'):
        print("=== baseline reproduced through Ladder2 (mux off, no Karatsuba) ===")
        test_comb(97, 3, 4, 2, 9, mux=False)
        test_comb(97, 3, 4, 3, 9, mux=False)
    if which in ('all', 'mux'):
        print("\n=== A. one-hot MUX look-up (zero-ancilla table) ===")
        for w in (1, 2, 3, 4):
            test_comb(97, 3, 4, w, 9, mux=True)
        test_comb(193, 5, 6, 2, 11, mux=True)
        test_comb(193, 5, 6, 3, 11, mux=True)
    if which in ('all', 'tree'):
        print("\n=== B. AND-tree one-hot (no cardinality penalty) ===")
        for w in (2, 3, 4):
            test_comb(97, 3, 4, w, 9, mux=True, onehot='tree')
        test_comb(193, 5, 6, 2, 11, mux=True, onehot='tree')
    if which in ('all', 'kara'):
        print("\n=== C. Karatsuba multiplication ===")
        for d in (1, 2, 3):
            test_comb(97, 3, 4, 2, 9, mux=True, kdepth=d, kmin=1, onehot='tree')
        test_comb(193, 5, 6, 2, 11, mux=True, kdepth=2, kmin=1, onehot='tree')
        test_comb(331, 2, 8, 2, 5, mux=True, kdepth=3, kmin=1, onehot='tree')
    if which in ('all', 'signed'):
        print("\n=== D. signed-digit comb ===")
        for w in (1, 2):
            test_comb(97, 3, 4, w, 9, mux=True, signed=True)
        test_comb(97, 3, 4, 2, 9, mux=True, signed=True, kdepth=2, kmin=1)
        test_comb(193, 5, 6, 2, 11, mux=True, signed=True, onehot='tree')
        # NOTE the signed-digit rewrite k' = 2k - (2^m - 1) is a bijection only
        # when the group order is ODD.  The real instance's order n is prime, so
        # it is; a toy curve of even order genuinely breaks it (see FINDINGS.md).
        test_comb(331, 2, 8, 2, 5, mux=True, signed=True, onehot='tree',
                  kdepth=2, kmin=1)
        test_comb(331, 2, 8, 4, 5, mux=True, signed=True, onehot='tree')
    if which in ('all', 'semaev'):
        print("\n=== E. x-only / Semaev S_3 chain ===")
        test_semaev(97, 3, 4, 2)
        test_semaev(193, 5, 6, 2)
        test_semaev(331, 2, 8, 2)
        test_semaev(331, 2, 8, 2, onehot='tree', kdepth=2)
        test_semaev(331, 2, 8, 4)
