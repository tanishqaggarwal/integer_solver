#!/usr/bin/env python3
"""crt.py -- residue-number-system / Freivalds arithmetization of one modular
multiplication, plus the exact cost of its base-conversion glue.

The relation to encode is always the same integer identity

        V  :=  A*B + p - W - p*Qw  ==  0            (A, B, W  s-bit words,
                                                     Qw an (s+1)-bit quotient)

Instead of balancing V column by column over Z (what qubo.assert_zero does),
check it modulo small coprime moduli m_1..m_r:

        A_i B_i + (p mod m_i) - W_i - (p mod m_i) Qw_i == 0   (mod m_i)

with X_i := X mod m_i.  |V| < 2^{2s+2} for any admissible words, so

    prod m_i > 2^{2s+3}   =>   the residue test is EXACT
    a single small modulus =>  Freivalds-style RELAXATION

The glue that makes either version work is BASE CONVERSION: X_i has to be
derived from X's BITS.  Nothing else bounds the size of an RNS-only value, and
without a size bound the residue identity is vacuous: for any A,B,W there is
always some Q with p*Q == A*B + p - W (mod M), since gcd(p, M) = 1.  So the
conversion cost is not optional overhead -- it is the whole scheme.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympy import nextprime               # noqa: E402
from enc import Ladder2                   # noqa: E402


def congruent_mod(L, poly, const, m, tag, valfn):
    """Ladder.congruent, but modulo m instead of p."""
    Q = L.qb
    lo = const + sum(min(0, c) for c in poly.values())
    hi = const + sum(max(0, c) for c in poly.values())
    qlo, qhi = lo // m, hi // m
    nb = max(0, (qhi - qlo).bit_length())
    Q.word(f"qm:{tag}", nb, lambda wv, f=valfn, b=qlo, m=m: f(wv) // m - b)
    poly = dict(poly)
    for t, v in enumerate(Q.trace[-1][2]):
        poly[(v,)] = poly.get((v,), 0) - m * (1 << t)
    Q.assert_zero(poly, const - m * qlo, tag)


def residue(L, bits, name, m, tag):
    """fresh word R == (value of `bits`) mod m, derived from the bits.
    THIS is the base conversion."""
    Q = L.qb
    R = Q.word(f"r:{tag}", m.bit_length(), lambda wv, n=name, m=m: wv[n] % m)
    poly = {}
    for t, v in enumerate(bits):
        poly[(v,)] = poly.get((v,), 0) + ((1 << t) % m)
    for t, v in enumerate(R):
        poly[(v,)] = poly.get((v,), 0) - (1 << t)
    nb = len(bits)
    congruent_mod(L, poly, 0, m, f"conv:{tag}",
                  lambda wv, n=name, m=m, nb=nb:
                  sum(((1 << t) % m) * ((wv[n] >> t) & 1) for t in range(nb))
                  - wv[n] % m)   # divisible by m by construction
    return R


def moduli_for(s, b, exact=True):
    """primes below 2^b, largest first, whose product exceeds 2^{2s+3} (exact),
    or the single largest one (Freivalds).

    Raises if the primes below 2^b cannot reach the bound -- which really can
    happen: there are only 54 primes below 2^8, and their product is ~2^369,
    short of the 2^515 an exact 256-bit check needs.  Small channels are not
    merely expensive, they are unavailable."""
    from sympy import primerange
    ps = sorted(primerange(2, 1 << b), reverse=True)
    if not exact:
        return ps[:1]
    out, prod, need = [], 1, 1 << (2 * s + 3)
    for q in ps:
        out.append(q)
        prod *= q
        if prod > need:
            return out
    raise ValueError(f"primes below 2^{b} multiply to only 2^{prod.bit_length()}"
                     f" < 2^{2*s+3}: an exact check is impossible at b={b}")


def build_rns_modmul(p, s, b, mode='binary', chunk=16, exact=True):
    """A*B == W (mod p) through r residue channels.  Witness inputs _A,_B,_W."""
    L = Ladder2(p, chunk=chunk, mode=mode)
    Q = L.qb
    A = Q.word("A", s, lambda wv: wv['_A'])
    Bw = Q.word("B", s, lambda wv: wv['_B'])
    W = Q.word("W", s, lambda wv: wv['_W'])
    Qw = Q.word("Qw", s + 1, lambda wv, p=p: wv['_Q'] if '_Q' in wv else
                (wv['_A'] * wv['_B'] + p - wv['_W']) // p)
    mods = moduli_for(s, b, exact)
    for i, m in enumerate(mods):
        Ai = residue(L, A, "A", m, f"A{i}")
        Bi = residue(L, Bw, "B", m, f"B{i}")
        Wi = residue(L, W, "W", m, f"W{i}")
        Qi = residue(L, Qw, "Qw", m, f"Q{i}")
        poly = {}
        for ii, u in enumerate(Ai):
            for jj, v in enumerate(Bi):
                mm = (u, v) if u <= v else (v, u)
                poly[mm] = poly.get(mm, 0) + (1 << (ii + jj))
        for t, v in enumerate(Wi):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        pm = p % m
        for t, v in enumerate(Qi):
            poly[(v,)] = poly.get((v,), 0) - pm * (1 << t)
        congruent_mod(
            L, poly, pm, m, f"chk{i}",
            lambda wv, m=m, pm=pm, p=p:
            (wv['_A'] % m) * (wv['_B'] % m) + pm - (wv['_W'] % m)
            - pm * ((wv['_Q'] if '_Q' in wv else
                     (wv['_A'] * wv['_B'] + p - wv['_W']) // p) % m))
    Q.finalize()
    return L, mods


def rns_modmul_cost(s, b, mode='binary', chunk=16, exact=True, seed=1):
    rnd = random.Random(seed)
    p = int(nextprime((1 << (s - 1)) + rnd.randrange(1 << (s - 2))))
    L, mods = build_rns_modmul(p, s, b, mode=mode, chunk=chunk, exact=exact)
    st = L.qb.stats()
    st['vars'] -= 3 * s          # A, B, W belong to the surrounding circuit
    st['nmods'] = len(mods)
    return st


def exhaustive(s, b, exact=True, p=None):
    """enumerate EVERY (A, B, W) triple of s-bit words AND every value of the
    free quotient word Qw, and compare the zero-energy set with
    { (A,B,W) : A*B == W (mod p) }.

    Qw is the one genuinely free ancilla here (every residue word is pinned by
    its conversion constraint), so letting it range is what makes this a real
    test of the relaxation rather than a test of the canonical witness."""
    p = p or int(nextprime(1 << (s - 1)))
    L, mods = build_rns_modmul(p, s, b, mode='wallace', exact=exact)
    Q = L.qb
    N = 1 << s
    zero = true = both = 0
    for a in range(N):
        for bb in range(N):
            for w in range(N):
                hit = False
                for q in range(1 << (s + 1)):
                    wv = {'_A': a, '_B': bb, '_W': w, '_Q': q}
                    try:
                        x, _ = Q.witness({}, wv)
                        e = Q.energy(x)
                    except Exception:
                        continue
                    if e == 0:
                        hit = True
                        break
                t = (a * bb - w) % p == 0
                zero += hit
                true += t
                both += (hit and t)
    return dict(p=p, mods=mods, vars=Q.n, zero=zero, true=true, both=both,
                spurious=zero - both, missed=true - both)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'cost'
    if what == 'check':
        print("exhaustive check of the CRT/Freivalds modmul gadget")
        for s, b, ex in ((4, 4, True), (4, 5, True), (4, 6, True),
                         (4, 3, False), (4, 4, False), (4, 5, False),
                         (4, 6, False), (4, 7, False)):
            r = exhaustive(s, b, ex)
            kind = "EXACT " if ex else "FREIVALDS"
            print(f"  s={s} b={b} {kind} r={len(r['mods'])} mods={r['mods']}: "
                  f"true={r['true']} zero-energy={r['zero']} "
                  f"spurious={r['spurious']} missed={r['missed']}")
    else:
        print("one 256x256 modular multiplication through RNS channels (exact)")
        print(f"{'b':>4} {'r':>4} {'vars':>10} {'couplers':>12} {'AND':>9} {'|J|':>6}")
        for b in (8, 12, 16, 24, 32, 48, 64):
            st = rns_modmul_cost(256, b)
            print(f"{b:4d} {st['nmods']:4d} {st['vars']:10,d} {st['couplers']:12,d} "
                  f"{st['and_vars']:9,d} 2^{st['dynamic_range_bits']}", flush=True)
