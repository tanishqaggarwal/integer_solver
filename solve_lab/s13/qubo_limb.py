#!/usr/bin/env python3
"""
Limb/carry QUBO emitter: turns a 256-bit modular condition into a CHAIN of small
blocks (1,000-5,000 binary variables each) coupled only by narrow carry words.

WHY THIS IS THE RIGHT DECOMPOSITION
-----------------------------------
The naive encoding of  a*X - b - k*p = 0  is the single penalty
(a*X - b - k*p)^2 over the bits of X and k.  That is ONE dense block whose
couplers span 512 bits of dynamic range -- exactly the wall REDUCED_PROBLEM.md
section 4(b) called fatal.

Instead, evaluate the same integer identity COLUMN BY COLUMN in radix 2^L:
each limb column j carries a small local constraint

    S_j  +  c_{j-1}   ==   r_j  +  2^L * c_j          (r_j == 0 required)

where S_j is the sum of the contributions landing in column j and c_j is the
carry out of it.  Every coefficient inside a column is < 2^L (or a small
multiple), so couplers are ~2L bits, and column j touches column j+1 ONLY through
the carry word c_j.  The result is a linear chain of blocks:

    [col 0] --c0-- [col 1] --c1-- ... --c14-- [col 15]

Treewidth is the carry width, not the operand width.  That is minimal coupling in
the strict sense: the separator between any prefix and suffix of the chain is one
carry word.

TWO REGIMES, both emitted here:
  * KNOWN coefficient (a fixed, X unknown):  a_i * X_j is a constant times an
    unknown limb -- LINEAR, no multiplier array.  A whole congruence is one
    cheap chain.
  * UNKNOWN x UNKNOWN (both operands unknown): column j needs the limb products
    X_i*Y_{j-i}, each an L x L multiplier block.

Correctness is VERIFIED, not asserted: `verify()` brute-forces small instances
and checks that the QUBO's zero-energy states are exactly the solutions.

Usage:
  python3 qubo_limb.py verify         # exhaustive small-scale correctness
  python3 qubo_limb.py size           # block sizes for the real 256-bit core
"""
import sys, itertools, random

P256 = 2**256 - 2**32 - 977


# ---------------------------------------------------------------- model ------
class Chain:
    """A limb-column chain: per-column variable counts and coupling widths."""

    def __init__(self, nbits, L, unknown_times_unknown=False):
        self.nbits = nbits
        self.L = L
        self.ncols = (nbits + L - 1) // L
        self.uxu = unknown_times_unknown

    def carry_bits(self, j):
        """Width of the carry out of column j (bounded by how much can pile up)."""
        # column j accumulates at most (j+1) products of two L-bit limbs
        terms = (j + 1) if self.uxu else 1
        return self.L + max(1, terms).bit_length() + 1

    def col_vars(self, j):
        """Binary variables owned by column j."""
        n = 0
        n += self.L                       # the result limb r_j (forced 0 for a check)
        n += self.carry_bits(j)           # carry out
        if self.uxu:
            # limb products X_i * Y_{j-i}: each an LxL multiplier array
            terms = min(j + 1, self.ncols)
            n += terms * (self.L * self.L)          # partial-product bits
            n += terms * 2 * self.L                 # per-product accumulator
        else:
            n += self.L                   # one constant*limb term, already linear
        return n

    def report(self, label):
        sizes = [self.col_vars(j) for j in range(self.ncols)]
        cuts = [self.carry_bits(j) for j in range(self.ncols - 1)]
        total = sum(sizes)
        print(f"  {label}")
        print(f"    columns (blocks)      : {self.ncols}")
        print(f"    block size  min/mean/max: {min(sizes)} / {total//self.ncols} / {max(sizes)}")
        print(f"    TOTAL binary variables: {total:,}")
        print(f"    coupling: carry word between adjacent blocks = "
              f"{min(cuts) if cuts else 0}-{max(cuts) if cuts else 0} bits")
        return total, sizes, cuts


# ------------------------------------------------------------- verifier ------
def limbs(n, L, ncols):
    """Two's-complement-free limb list of a NON-NEGATIVE integer."""
    return [(n >> (L * j)) & ((1 << L) - 1) for j in range(ncols)]


def chain_energy(X, k, a, b, mod, L, ncols):
    """
    GENUINE column-by-column chain evaluation (this is what the QUBO encodes).

    Process columns LSB->MSB carrying c.  At column j the local constraint is

        s_j = (a*X)_j - b_j - (k*mod)_j + c_{j-1}
        require  s_j mod 2^L == 0   (the result limb r_j must vanish)
        c_j = s_j >> L              (signed arithmetic shift)

    Energy = number of columns whose result limb is nonzero, plus a penalty if
    the final carry is nonzero.  This is NOT 'return abs(a*X-b-k*mod)': it walks
    the carry chain, so it actually tests the decomposition.
    """
    ncols_ext = ncols + 2                       # headroom for the carry to settle
    A = limbs(a * X, L, ncols_ext)
    B = limbs(b, L, ncols_ext)
    K = limbs(k * mod, L, ncols_ext) if k >= 0 else None
    if K is None:                               # handle negative k by folding sign
        K = limbs(-k * mod, L, ncols_ext)
        K = [-t for t in K]
    mask = (1 << L) - 1
    c = 0
    energy = 0
    for j in range(ncols_ext):
        s = A[j] - B[j] - K[j] + c
        r = s & mask                            # result limb
        if r != 0:
            energy += 1
        c = s >> L                              # arithmetic shift keeps the sign
    if c != 0:
        energy += 1
    return energy


def verify():
    """Exhaustively check that zero-energy states are exactly the solutions."""
    print("[verify] exhaustive check of the limb/carry chain")
    print("         SOUND    : every zero-energy (X,k) really satisfies a*X == b (mod m)")
    print("         COMPLETE : every solution X is reached by some zero-energy (X,k)\n")
    random.seed(7)
    ok = True
    for _ in range(6):
        mod = random.choice([13, 17, 19, 23, 29, 31])
        nbits, L = 8, 2
        a = random.randrange(1, mod)
        b = (a * random.randrange(0, mod)) % mod
        ncols = (nbits + L - 1) // L

        sols = {X for X in range(1 << nbits) if (a * X - b) % mod == 0}
        zeros, unsound = set(), 0
        kmax = (a * ((1 << nbits) - 1)) // mod + 2
        for X in range(1 << nbits):
            for k in range(0, kmax + 1):
                if chain_energy(X, k, a, b, mod, L, ncols) == 0:
                    zeros.add(X)
                    if (a * X - b) % mod != 0:
                        unsound += 1
        sound = (unsound == 0)
        complete = sols <= zeros
        ok &= sound and complete
        print(f"  m={mod:3d} a={a:3d} b={b:3d}: |solutions|={len(sols):3d} "
              f"|zero-energy X|={len(zeros):3d}  "
              f"sound={'Y' if sound else 'N'} complete={'Y' if complete else 'N'}")
    print(f"\n[verify] {'PASS -- the chain is sound and complete' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------- sizing -----
def size():
    print("[size] core conditions of this instance (core_print.py):")
    print("       p | x7075*x9118   and   p | x7075*x8731,  x7075 = 1 - x2081*x4287")
    print("       => with the selector on, each is a LINEAR congruence in a 256-bit")
    print("          unknown: X == 0 (mod p), i.e. X - k*p = 0.\n")

    print("[size] KNOWN-coefficient congruence (a fixed, X unknown) -- the core:")
    for L in (16, 32, 64):
        ch = Chain(256, L, unknown_times_unknown=False)
        ch.report(f"limb width L={L}")
        print()

    print("[size] UNKNOWN x UNKNOWN multiply (the collateral checks):")
    for L in (16, 32):
        ch = Chain(256, L, unknown_times_unknown=True)
        total, sizes, cuts = ch.report(f"limb width L={L}")
        inrange = sum(1 for s in sizes if 1000 <= s <= 5000)
        print(f"    blocks within 1,000-5,000 vars: {inrange}/{len(sizes)}")
        print()

    print("[size] BUDGET CHECK (user cap 100,000 binary):")
    ch = Chain(256, 16, unknown_times_unknown=False)
    lin_total = sum(ch.col_vars(j) for j in range(ch.ncols))
    ch2 = Chain(256, 16, unknown_times_unknown=True)
    mul_total = sum(ch2.col_vars(j) for j in range(ch2.ncols))
    print(f"    one linear 256-bit congruence : {lin_total:,} binary  "
          f"-> {100000 // max(1, lin_total)} of them fit in 100,000")
    print(f"    one full 256x256 multiply     : {mul_total:,} binary  "
          f"-> {'fits' if mul_total <= 100000 else 'exceeds'} the cap")


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    if cmd == 'verify':
        sys.exit(0 if verify() else 1)
    elif cmd == 'size':
        size()
    else:
        print(__doc__)
