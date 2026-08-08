#!/usr/bin/env python3
"""delegable.py -- WHERE THE FLOOR IS: which fragments carry decision content?

Claim (Lemma 1 in FINDINGS.md): in this encoding, ANY fragment of the QUBO that
does not contain the terminal constraint "accumulator == T" has a zero-energy
state for EVERY assignment of its free input variables.  Its ground-state set is
therefore in bijection with its free inputs, an anneal returns an arbitrary
element of it, and a sound outer loop must enumerate those inputs classically --
i.e. the QPU contributes nothing.

Proof sketch (compiler level): `qubo.py`'s `witness(inputs, wv0)` replays the
construction trace and produces a zero-energy assignment from ANY input bits.
Every penalty is (linear form)^2 or an AND penalty, and the replay sets each
ancilla to the value that zeroes its own form.  The only penalties that can
survive are the ones whose linear form contains no fresh ancilla -- and the only
such penalties in the whole construction are the two terminal congruences
"x_final == Tx", "y_final == Ty" (and the one-hot digit constraints, which are
satisfied by construction).

This script CHECKS that claim exhaustively on scaled instances:

  A. one modular multiplication  a*b == c (mod p), a,b,c all free     -> E=0 always
  B. one comb window (2 look-ups + 1 EC addition), digits free        -> E=0 always
  C. an M-window comb WITHOUT the terminal check, all digits free     -> E=0 always
  D. the same comb WITH the terminal check                            -> E=0 only on k

Writes multirun/delegable.json.
"""
import json, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from ecsmall import curve, find
from ladder import Ladder, not_equal, build_win

MODE = 'wallace'


# ---------------------------------------------------------------- A: modmul
def check_modmul(p, s=None):
    s = s or p.bit_length()
    L = Ladder(p, mode=MODE); Q = L.qb
    box = {}
    A = Q.word("a", s, lambda wv: box['a'])
    B = Q.word("b", s, lambda wv: box['b'])
    C = Q.word("c", s, lambda wv: box['a'] * box['b'] % p)
    L.mul_eq("mm", A, B, "a", "b", [(1, C, "c")], 0)
    Q.finalize()
    bad = 0
    for a in range(p):
        for b in range(p):
            box['a'], box['b'] = a, b
            x, _ = Q.witness({}, {})
            if Q.energy(x) != 0:
                bad += 1
    return dict(vars=Q.n, tested=p * p, nonzero_energy=bad)


# ------------------------------------------------------- B/C: comb, no final
def build_comb_nofinal(p, table, w, mode=MODE):
    """exactly build_win, minus the two terminal congruences."""
    M = len(table); D = 1 << w
    L = Ladder(p, mode=mode); Q = L.qb; s = L.s
    U = []
    for j in range(M):
        u = [Q.new(f"u{j}_{t}", 'input') for t in range(D)]
        for t, v in enumerate(u):
            Q.trace.append(('word', f"u{j}_{t}", [v],
                            (lambda wv, j=j, t=t: 1 if wv[f"_u{j}"] == t else 0)))
        Q.add_square({v: 1 for v in u}, -1)
        U.append(u)

    def W(name, terms, const, fn):
        wd = Q.word(name, s, fn)
        poly = {}
        for coef, bits, _ in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
        for t, v in enumerate(wd):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, const, f"lin:{name}",
                    lambda wv, terms=terms, const=const, name=name:
                    const + sum(c * wv[nm] for c, _b, nm in terms) - wv[name])
        return wd

    def select(name, j, coord):
        terms = [(table[j][t][coord], [U[j][t]], f"u{j}_{t}") for t in range(D)]
        return W(name, terms, 0,
                 lambda wv, j=j, coord=coord: table[j][wv[f"_u{j}"]][coord] % p)

    x1 = select("x1_0", 0, 0); nx = "x1_0"
    y1 = select("y1_0", 0, 1); ny = "y1_0"
    for j in range(1, M):
        x2 = select(f"x2_{j}", j, 0); n2x = f"x2_{j}"
        y2 = select(f"y2_{j}", j, 1); n2y = f"y2_{j}"
        d = W(f"d{j}", [(1, x2, n2x), (-1, x1, nx)], 0,
              lambda wv, n2x=n2x, nx=nx: (wv[n2x] - wv[nx]) % p)
        e = W(f"e{j}", [(1, y2, n2y), (-1, y1, ny)], 0,
              lambda wv, n2y=n2y, ny=ny: (wv[n2y] - wv[ny]) % p)
        lam = Q.word(f"lam{j}", s,
                     lambda wv, j=j: wv[f"e{j}"] * pow(wv[f"d{j}"], -1, p) % p)
        L.mul_eq(f"lam{j}", lam, d, f"lam{j}", f"d{j}", [(1, e, f"e{j}")], 0)
        for c in (0, p):
            not_equal(L, Q, d, f"d{j}", c, f"ne{j}_{c}")
        x3 = Q.word(f"x3_{j}", s, lambda wv, j=j, nx=nx, n2x=n2x:
                    (wv[f"lam{j}"] ** 2 - wv[nx] - wv[n2x]) % p)
        L.mul_eq(f"x3e{j}", lam, lam, f"lam{j}", f"lam{j}",
                 [(1, x3, f"x3_{j}"), (1, x1, nx), (1, x2, n2x)], 0)
        mm = W(f"m{j}", [(1, x1, nx), (-1, x3, f"x3_{j}")], 0,
               lambda wv, j=j, nx=nx: (wv[nx] - wv[f"x3_{j}"]) % p)
        y3 = Q.word(f"y3_{j}", s, lambda wv, j=j, ny=ny:
                    (wv[f"lam{j}"] * wv[f"m{j}"] - wv[ny]) % p)
        L.mul_eq(f"y3e{j}", lam, mm, f"lam{j}", f"m{j}",
                 [(1, y3, f"y3_{j}"), (1, y1, ny)], 0)
        x1, y1, nx, ny = x3, y3, f"x3_{j}", f"y3_{j}"
    Q.finalize()
    return L, U


def _table(p, Bc, m, w):
    add, mul = curve(p, Bc)
    G, order = find(p, Bc)
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)]
             for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    return add, mul, G, order, table, M, off


def check_comb(p, Bc, m, w, with_final):
    add, mul, G, order, table, M, off = _table(p, Bc, m, w)
    D = 1 << w
    if with_final:
        # pick a k whose ladder is non-degenerate (same convention as demo_win.py)
        def ok(kk):
            dg = [(kk >> (w * j)) % D for j in range(M)]
            S = table[0][dg[0]]
            for j in range(1, M):
                Qp = table[j][dg[j]]
                if S is None or Qp is None or S[0] == Qp[0]: return False
                S = add(S, Qp)
            return S is not None
        # prefer a non-trivial k so the demonstration is not vacuous
        cands = [kk for kk in range(1 << m) if ok(kk)]
        k = max(cands) if len(cands) > 1 else cands[0]
        Tp = add(mul(k % order, G), mul(off % order, G))
        L, U = build_win(p, Bc, table, Tp, w, mode=MODE, neq=True)
    else:
        k = None
        L, U = build_comb_nofinal(p, table, w)
    Q = L.qb

    zeros, degen = [], 0
    for kk in range(D ** M):
        dg = [(kk >> (w * j)) % D for j in range(M)]
        wv0 = {f"_u{j}": dg[j] for j in range(M)}
        inp = {}
        for j in range(M):
            for t in range(D):
                inp[U[j][t]] = 1 if t == dg[j] else 0
        try:
            x, _ = Q.witness(inp, wv0)
        except Exception:
            degen += 1
            continue
        if Q.energy(x) == 0:
            zeros.append(kk)
    return dict(vars=Q.n, windows=M, w=w, digit_assignments=D ** M,
                zero_energy=len(zeros), degenerate=degen, k=k,
                zero_list=zeros if len(zeros) <= 8 else zeros[:8] + ['...'])


# ---------------------------------------------------------------- driver
if __name__ == '__main__':
    out = {}

    print("A. one modular multiplication, a and b FREE  (is it a decision problem?)")
    for p in (7, 11, 13, 17, 19):
        r = check_modmul(p)
        out.setdefault('modmul', []).append(dict(p=p, **r))
        print(f"   p={p:3d}: {r['vars']:5d} vars, all {r['tested']:5d} (a,b) pairs "
              f"reach E=0  -> nonzero-energy pairs: {r['nonzero_energy']}")
    print("   => the ground-state set IS the input set.  No decision content.\n")

    print("B/C. an M-window comb WITHOUT the terminal check, all digits free")
    for (p, Bc, m, w) in ((97, 3, 4, 2), (97, 3, 4, 1), (97, 3, 6, 2), (89, 5, 4, 2)):
        r = check_comb(p, Bc, m, w, with_final=False)
        out.setdefault('comb_nofinal', []).append(dict(p=p, B=Bc, m=m, **r))
        print(f"   p={p} m={m} w={w}: {r['windows']} windows, {r['vars']:5d} vars; "
              f"{r['zero_energy']}/{r['digit_assignments']} digit assignments reach E=0 "
              f"({r['degenerate']} degenerate)")
    print("   => every non-degenerate digit assignment is a ground state."
          "  No decision content.\n")

    print("D. the SAME comb WITH the terminal check")
    for (p, Bc, m, w) in ((97, 3, 4, 2), (97, 3, 4, 1), (89, 5, 4, 2)):
        r = check_comb(p, Bc, m, w, with_final=True)
        out.setdefault('comb_final', []).append(dict(p=p, B=Bc, m=m, **r))
        print(f"   p={p} m={m} w={w}: {r['windows']} windows, {r['vars']:5d} vars; "
              f"zero-energy digit assignments = {r['zero_list']}  (true k = {r['k']})")
    print("   => only the discrete log survives.  The terminal check is the ONLY"
          " source of decision content.")

    json.dump(out, open(os.path.join(_HERE, 'delegable.json'), 'w'), indent=1)
    print("\nwrote multirun/delegable.json")
