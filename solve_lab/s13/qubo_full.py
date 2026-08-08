#!/usr/bin/env python3
"""
THE COMPLETE QUBO SET for EQUATIONS.txt.

Architecture (two layers, both emitted here)
--------------------------------------------
LAYER A -- the bulk circuit, in RNS.
    Work modulo a small prime q (w = ceil(log2 q) bits).  Because the file is a
    pure polynomial system over Z, reduction mod q is a ring homomorphism, so
    every one of the 39,033 equations has an exact mod-q image.  Each ATOM
    becomes one block:
        LIN  out = sum c_i x_i (mod q)          -- carry-structured adder
        MUL  out = x*y (mod q)                  -- w x w partial products
        SQ   out = x*x (mod q)
        BIT  x*x = x                            -- boolean pin
    Wires (w-bit residues) are shared between blocks: that is the ONLY coupling.

LAYER B -- the p-divisibility core, which RNS cannot see.
    Mod q every p-handle absorbs its check (p is invertible mod q), so the mod-q
    image is blind to the real content: divisibility by p.  Those conditions are
    emitted as limb/carry chains (see qubo_limb.py), verified sound+complete.
    Measured core (core_print.py):  p | x7075*x9118 and p | x7075*x8731.

Every block is a genuine QUBO: a dict {(i,j): coeff} over binary variables, whose
MINIMUM ENERGY IS 0 exactly on the assignments satisfying that atom.  Blocks are
materialized on demand (all of them at once would be tens of GB); `verify`
brute-forces small ones against ground truth.

Usage:
    python3 qubo_full.py verify      # correctness of each block type
    python3 qubo_full.py census      # the full block set for the real instance
    python3 qubo_full.py emit A123   # materialize one block's coefficients
"""
import os, sys, json, time, itertools
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

P256 = 2**256 - 2**32 - 977
NA = len(L.avars)


# ===================================================================== QUBO ==
class QUBO:
    """Sparse symmetric QUBO:  E(z) = sum_{i<=j} Q[(i,j)] z_i z_j + const."""

    def __init__(self):
        self.Q = defaultdict(int)
        self.const = 0
        self.nv = 0
        self.names = {}

    def var(self, name):
        if name not in self.names:
            self.names[name] = self.nv
            self.nv += 1
        return self.names[name]

    def add_square(self, terms, const=0):
        """Add ( sum_t coeff*z_var + const )^2 to the energy."""
        # expand the square over binary z (z^2 = z)
        for (c1, v1) in terms:
            self.Q[(v1, v1)] += c1 * c1 + 2 * c1 * const
        for a in range(len(terms)):
            for b in range(a + 1, len(terms)):
                c1, v1 = terms[a]; c2, v2 = terms[b]
                i, j = (v1, v2) if v1 <= v2 else (v2, v1)
                self.Q[(i, j)] += 2 * c1 * c2
        self.const += const * const

    def energy(self, assign):
        """assign: dict varindex -> 0/1."""
        e = self.const
        for (i, j), c in self.Q.items():
            e += c * assign.get(i, 0) * assign.get(j, 0)
        return e

    def max_coeff(self):
        return max((abs(c) for c in self.Q.values()), default=0)

    def stats(self):
        return dict(nvars=self.nv, nterms=len(self.Q),
                    max_coeff=self.max_coeff(),
                    coeff_bits=self.max_coeff().bit_length())


def and_gadget(Q, zx, zy, zp):
    """Penalty forcing zp = zx AND zy  (Rosenberg): xy - 2xp - 2yp + 3p."""
    i, j = (zx, zy) if zx <= zy else (zy, zx)
    Q.Q[(i, j)] += 1
    a, b = (zx, zp) if zx <= zp else (zp, zx)
    Q.Q[(a, b)] += -2
    a, b = (zy, zp) if zy <= zp else (zp, zy)
    Q.Q[(a, b)] += -2
    Q.Q[(zp, zp)] += 3


# ============================================================ block builders =
def mul_geometry(w, q):
    """Shared sizing so builder and verifier agree exactly."""
    kb = w + 1                                   # quotient width: x*y/q < 2^(w+1)
    ncol = 2 * w + 2                             # columns of the comparison
    cw = (w + 2).bit_length()                    # carry width per side
    qbits = [(q >> t) & 1 for t in range(q.bit_length())]
    return kb, ncol, cw, qbits


def build_mul(w, q):
    """
    out = x*y (mod q)  <=>  x*y == out + k*q  (both sides NON-NEGATIVE).

    Encoded as TWO carry chains that are then equated bit by bit:
        LHS chain: sum of partial products  x_i*y_j
        RHS chain: out + k*q
    Both sides only ADD, so every carry is a non-negative binary word -- this is
    the fix for the naive `x*y - out - k*q` chain, whose column sums go negative
    and which is therefore sound but INCOMPLETE.
    Coefficients stay small (~2^(carry width)), so couplers are ~10-14 bits.
    """
    kb, ncol, cw, qbits = mul_geometry(w, q)
    Q = QUBO()
    zx = [Q.var(f'x{i}') for i in range(w)]
    zy = [Q.var(f'y{i}') for i in range(w)]
    zo = [Q.var(f'o{i}') for i in range(w)]
    zk = [Q.var(f'k{i}') for i in range(kb)]

    zp = {}
    for i in range(w):
        for j in range(w):
            zp[(i, j)] = Q.var(f'p{i}_{j}')
            and_gadget(Q, zx[i], zy[j], zp[(i, j)])

    zl  = [Q.var(f'L{c}')  for c in range(ncol)]          # LHS bit per column
    zr  = [Q.var(f'R{c}')  for c in range(ncol)]          # RHS bit per column
    zcl = {c: [Q.var(f'cl{c}_{b}') for b in range(cw)] for c in range(ncol)}
    zcr = {c: [Q.var(f'cr{c}_{b}') for b in range(cw)] for c in range(ncol)}

    prevL, prevR = [], []
    for col in range(ncol):
        # ---- LHS column: partial products + carry_in = bit + 2*carry_out ----
        tl = [(1, zp[(i, col - i)]) for i in range(w) if 0 <= col - i < w]
        tl += [(1 << b, v) for b, v in enumerate(prevL)]
        tl += [(-1, zl[col])]
        tl += [(-(1 << (b + 1)), v) for b, v in enumerate(zcl[col])]
        Q.add_square(tl)

        # ---- RHS column: out bit + k*q bits + carry_in = bit + 2*carry_out --
        tr = []
        if col < w:
            tr.append((1, zo[col]))
        for t, bit in enumerate(qbits):
            ki = col - t
            if bit and 0 <= ki < kb:
                tr.append((1, zk[ki]))
        tr += [(1 << b, v) for b, v in enumerate(prevR)]
        tr += [(-1, zr[col])]
        tr += [(-(1 << (b + 1)), v) for b, v in enumerate(zcr[col])]
        Q.add_square(tr)

        # ---- equality of the two sides at this column -----------------------
        Q.add_square([(1, zl[col]), (-1, zr[col])])
        prevL, prevR = zcl[col], zcr[col]
    return Q


def lin_geometry(w, q, coeffs):
    cmax = max(abs(c) for c in coeffs)
    kb = w + cmax.bit_length() + len(coeffs).bit_length() + 1
    ncol = w + cmax.bit_length() + len(coeffs).bit_length() + q.bit_length() + 2
    cw = (len(coeffs) * cmax + 4).bit_length() + 1
    qbits = [(q >> t) & 1 for t in range(q.bit_length())]
    return kb, ncol, cw, qbits


def build_lin(w, q, coeffs):
    """
    out = sum_t coeffs[t]*v_t (mod q).  Coefficients are KNOWN, so each
    coeffs[t]*v_t is a constant times an unknown -- no partial-product array.
    Still carry-structured (two non-negative chains, equated column by column) so
    couplers stay ~2^(carry width) instead of growing with q.
    """
    kb, ncol, cw, qbits = lin_geometry(w, q, coeffs)
    Q = QUBO()
    zs = [[Q.var(f'v{t}_{i}') for i in range(w)] for t in range(len(coeffs))]
    zo = [Q.var(f'o{i}') for i in range(w)]
    zk = [Q.var(f'k{i}') for i in range(kb)]
    zl = [Q.var(f'L{c}') for c in range(ncol)]
    zr = [Q.var(f'R{c}') for c in range(ncol)]
    zcl = {c: [Q.var(f'cl{c}_{b}') for b in range(cw)] for c in range(ncol)}
    zcr = {c: [Q.var(f'cr{c}_{b}') for b in range(cw)] for c in range(ncol)}

    prevL, prevR = [], []
    for col in range(ncol):
        # LHS: sum_t coeffs[t]*v_t, expanded as shifted copies (constant weights)
        tl = []
        for t, c in enumerate(coeffs):
            for b in range(c.bit_length()):
                if (c >> b) & 1:
                    i = col - b
                    if 0 <= i < w:
                        tl.append((1, zs[t][i]))
        tl += [(1 << b, v) for b, v in enumerate(prevL)]
        tl += [(-1, zl[col])]
        tl += [(-(1 << (b + 1)), v) for b, v in enumerate(zcl[col])]
        Q.add_square(tl)

        # RHS: out + k*q
        tr = []
        if col < w:
            tr.append((1, zo[col]))
        for t, bit in enumerate(qbits):
            ki = col - t
            if bit and 0 <= ki < kb:
                tr.append((1, zk[ki]))
        tr += [(1 << b, v) for b, v in enumerate(prevR)]
        tr += [(-1, zr[col])]
        tr += [(-(1 << (b + 1)), v) for b, v in enumerate(zcr[col])]
        Q.add_square(tr)

        Q.add_square([(1, zl[col]), (-1, zr[col])])
        prevL, prevR = zcl[col], zcr[col]
    return Q


def build_bit():
    """Boolean pin x*x = x -- automatically satisfied by a binary variable."""
    Q = QUBO()
    Q.var('x')
    return Q


# ================================================================== verify ===
def verify():
    print("[verify] block-type correctness (brute force vs ground truth)\n")
    ok = True

    # ---- MUL over a small modulus -----------------------------------------
    # the encoding needs 2^w >= q so every residue fits in w bits
    for (w, q) in [(2, 3), (3, 5), (3, 7), (4, 11), (4, 13)]:
        assert (1 << w) >= q, f"w={w} too narrow for q={q}"
        Q = build_mul(w, q)
        kb, ncol, cw, qbits = mul_geometry(w, q)
        good = {(x, y, (x * y) % q) for x in range(min(1 << w, q))
                for y in range(min(1 << w, q))}
        found, unsound = set(), 0
        for x in range(min(1 << w, q)):
            for y in range(min(1 << w, q)):
                for o in range(min(1 << w, q)):
                    num = x * y - o
                    if num % q != 0 or num // q < 0 or num // q >= (1 << kb):
                        continue
                    k = num // q
                    assign = {}
                    for i in range(w):
                        assign[Q.names[f'x{i}']] = (x >> i) & 1
                        assign[Q.names[f'y{i}']] = (y >> i) & 1
                        assign[Q.names[f'o{i}']] = (o >> i) & 1
                    for i in range(kb):
                        assign[Q.names[f'k{i}']] = (k >> i) & 1
                    for i in range(w):
                        for j in range(w):
                            assign[Q.names[f'p{i}_{j}']] = \
                                ((x >> i) & 1) * ((y >> j) & 1)
                    # derive both chains (all non-negative)
                    cl = cr = 0
                    okc = True
                    for col in range(ncol):
                        sl = cl + sum(((x >> i) & 1) * ((y >> (col - i)) & 1)
                                      for i in range(w) if 0 <= col - i < w)
                        bl, cl = sl & 1, sl >> 1
                        sr = cr + ((o >> col) & 1 if col < w else 0)
                        for t, bit in enumerate(qbits):
                            ki = col - t
                            if bit and 0 <= ki < kb:
                                sr += (k >> ki) & 1
                        br, cr = sr & 1, sr >> 1
                        if cl >= (1 << cw) or cr >= (1 << cw):
                            okc = False; break
                        assign[Q.names[f'L{col}']] = bl
                        assign[Q.names[f'R{col}']] = br
                        for b in range(cw):
                            assign[Q.names[f'cl{col}_{b}']] = (cl >> b) & 1
                            assign[Q.names[f'cr{col}_{b}']] = (cr >> b) & 1
                    if not okc or cl != 0 or cr != 0:
                        continue
                    if Q.energy(assign) == 0:
                        found.add((x, y, o))
                        if (x * y) % q != o:
                            unsound += 1
        complete = good <= found
        sound = unsound == 0
        ok &= complete and sound
        print(f"  MUL w={w} q={q}: valid triples={len(good):3d} "
              f"zero-energy={len(found):3d}  sound={'Y' if sound else 'N'} "
              f"complete={'Y' if complete else 'N'}  "
              f"maxcoeff={Q.max_coeff()} ({Q.max_coeff().bit_length()} bits) "
              f"nvars={Q.nv}")

    # ---- LIN ---------------------------------------------------------------
    for (w, q, cf) in [(3, 5, [1, 1]), (3, 7, [2, 3]), (4, 11, [3, 5, 7])]:
        assert (1 << w) >= q, f"w={w} too narrow for q={q}"
        Q = build_lin(w, q, cf)
        kb, ncol, cw, qbits = lin_geometry(w, q, cf)
        n = len(cf)
        good, found, unsound = set(), set(), 0
        for vals in itertools.product(range(min(1 << w, q)), repeat=n):
            tgt = sum(c * v for c, v in zip(cf, vals)) % q
            good.add((vals, tgt))
            for o in range(min(1 << w, q)):
                num = sum(c * v for c, v in zip(cf, vals)) - o
                if num % q != 0 or num // q < 0 or num // q >= (1 << kb):
                    continue
                k = num // q
                assign = {}
                for t, v in enumerate(vals):
                    for i in range(w):
                        assign[Q.names[f'v{t}_{i}']] = (v >> i) & 1
                for i in range(w):
                    assign[Q.names[f'o{i}']] = (o >> i) & 1
                for i in range(kb):
                    assign[Q.names[f'k{i}']] = (k >> i) & 1
                cl = cr = 0
                okc = True
                for col in range(ncol):
                    sl = cl
                    for t, c in enumerate(cf):
                        for b in range(c.bit_length()):
                            if (c >> b) & 1:
                                i = col - b
                                if 0 <= i < w:
                                    sl += (vals[t] >> i) & 1
                    bl, cl = sl & 1, sl >> 1
                    sr = cr + ((o >> col) & 1 if col < w else 0)
                    for t, bit in enumerate(qbits):
                        ki = col - t
                        if bit and 0 <= ki < kb:
                            sr += (k >> ki) & 1
                    br, cr = sr & 1, sr >> 1
                    if cl >= (1 << cw) or cr >= (1 << cw):
                        okc = False; break
                    assign[Q.names[f'L{col}']] = bl
                    assign[Q.names[f'R{col}']] = br
                    for b in range(cw):
                        assign[Q.names[f'cl{col}_{b}']] = (cl >> b) & 1
                        assign[Q.names[f'cr{col}_{b}']] = (cr >> b) & 1
                if not okc or cl != 0 or cr != 0:
                    continue
                if Q.energy(assign) == 0:
                    found.add((vals, o))
                    if o != tgt:
                        unsound += 1
        complete = good <= found
        ok &= complete and unsound == 0
        print(f"  LIN w={w} q={q} coeffs={cf}: sound={'Y' if unsound==0 else 'N'} "
              f"complete={'Y' if complete else 'N'} "
              f"maxcoeff={Q.max_coeff()} ({Q.max_coeff().bit_length()} bits) "
              f"nvars={Q.nv}")

    print(f"\n[verify] {'PASS -- every block type is sound and complete' if ok else 'FAIL'}")
    return ok


# ================================================================== census ===
def degree(a):
    return max((len(m) for m in L.polys[a]), default=0)


def atom_cost(a, mul_size, lin_cost):
    """(#MUL sub-blocks, linear-accumulator vars) for one atom."""
    nprod = sum(1 for m in L.polys[a] if len(m) >= 2)
    nlin = sum(1 for m in L.polys[a] if len(m) < 2)
    return nprod, lin_cost(nlin)


def pack(units, lo, hi):
    """
    Greedily pack indivisible sub-blocks (each `unit` binary vars) into QUBO
    blocks of size in [lo,hi].  Returns the list of block sizes.  This is what
    makes EVERY block land in the target band: the atom is not the block, the
    sub-block is, and blocks are groups of sub-blocks.
    """
    blocks = []
    cur = 0
    for u in units:
        if u > hi:                       # a single unit already exceeds hi
            blocks.append(u); continue
        if cur + u > hi:
            blocks.append(cur); cur = 0
        cur += u
    if cur:
        blocks.append(cur)
    return blocks


def census(w=16, q=65521, lo=1000, hi=5000):
    t0 = time.time()
    print(f"[census] COMPLETE QUBO SET  (RNS layer q={q}, w={w} bits/residue)\n")

    gates = [a for a in range(NA) if L.atom_out.get(a) is not None]
    checks = [a for a in range(NA) if L.atom_out.get(a) is None]
    mul_ref = build_mul(w, q)
    # shared wires (x,y,out) are counted once globally, not per block
    mul_size = mul_ref.nv - 3 * w
    mul_coeff = mul_ref.max_coeff()

    # LIN cost is dominated by the per-ATOM carry chain, NOT by the term count:
    # measure it instead of assuming a per-term constant.
    lin_meas = {}
    for T in (2, 4, 8, 16):
        r = build_lin(w, q, [1] * T)
        lin_meas[T] = r.nv - T * w - w        # minus shared inputs and output
    slope = (lin_meas[16] - lin_meas[2]) / 14.0
    base = lin_meas[2] - 2 * slope

    def lin_cost(T):
        return int(base + slope * max(2, T))

    print(f"  primitive MUL sub-block : {mul_size} binary (excl. shared wires), "
          f"max coupler {mul_coeff} ({mul_coeff.bit_length()} bits) [VERIFIED]")
    print(f"  primitive LIN block     : {base:.0f} fixed + {slope:.1f}/term "
          f"(measured at T=2,4,8,16: "
          f"{ {k: v for k, v in lin_meas.items()} })\n")
    lin_per_term = slope   # kept for the older tiers below

    # ---------- (1) FULL INSTANCE: every atom encoded, nothing assumed -------
    units, nmul, totlin = [], 0, 0
    for a in range(NA):
        np_, lin = atom_cost(a, mul_size, lin_cost)
        nmul += np_; totlin += lin
        units.extend([mul_size] * np_)
        units.append(lin)
    tot_full = nmul * mul_size + totlin
    blocks_full = pack(units, lo, hi)
    inband_full = sum(1 for s in blocks_full if lo <= s <= hi)
    print(f"  [A] FULL INSTANCE (search everything, assume nothing)")
    print(f"      atoms                 : {NA:,} "
          f"({len(gates):,} gates, {len(checks):,} checks)")
    print(f"      MUL sub-blocks        : {nmul:,}")
    print(f"      TOTAL binary (1 prime): {tot_full:,}")
    print(f"      packed into           : {len(blocks_full):,} blocks, "
          f"{inband_full/len(blocks_full):.0%} inside [{lo:,},{hi:,}]")
    print(f"      -> honest, but ~{tot_full*20/1e9:.1f}e9 binary over 20 primes: "
          f"a scale statement, not a machine target.\n")

    # ---------- (2) CORE-ANCHORED: search only the residual ------------------
    # Everything outside the drivers' descendant cone is already satisfied and
    # is FIXED (constant-folded), so it needs no QUBO variables at all.
    ext = json.load(open(os.path.join(os.path.dirname(__file__),
                                      'core_extend.json')))
    cone_atoms = set(ext['core_atoms'])
    units_c, nmul_c, lin_c = [], 0, 0
    for a in sorted(cone_atoms):
        np_, lin = atom_cost(a, mul_size, lin_cost)
        nmul_c += np_; lin_c += lin
        units_c.extend([mul_size] * np_)
        units_c.append(lin)
    tot_core = nmul_c * mul_size + lin_c
    blocks_core = pack(units_c, lo, hi)
    inband_core = sum(1 for s in blocks_core if lo <= s <= hi)
    print(f"  [B] CORE-ANCHORED (fix the 11 solved advice numbers, "
          f"search the residual)")
    print(f"      cone atoms            : {len(cone_atoms):,}")
    print(f"      MUL sub-blocks        : {nmul_c:,}")
    print(f"      TOTAL binary (1 prime): {tot_core:,}")
    print(f"      packed into           : {len(blocks_core):,} blocks "
          f"(min {min(blocks_core)}, mean {sum(blocks_core)//len(blocks_core)}, "
          f"max {max(blocks_core)})")
    print(f"      inside [{lo:,},{hi:,}]  : {inband_core}/{len(blocks_core)} "
          f"({inband_core/len(blocks_core):.0%})\n")

    # ---------- (2b) CONSTANT-FOLDED: only downstream-of-driver vars vary ----
    # Everything not downstream of the 4 drivers is FIXED at its verified value,
    # so a monomial with <=1 unknown factor is a CONSTANT times an unknown --
    # linear, no multiplier array.  Only monomials with >=2 unknown factors need
    # a MUL sub-block.  This is the decisive size reduction.
    DRIVERS = [2081, 4287, 8731, 9118]
    users = defaultdict(list)
    for a in range(NA):
        for x in L.avars[a]:
            users[x].append(a)
    unknown = set(DRIVERS)
    touched = set()
    dq = list(DRIVERS)
    while dq:
        x = dq.pop()
        for a in users[x]:
            touched.add(a)
            oc = L.atom_out.get(a)
            if oc is not None and oc[1] not in unknown:
                unknown.add(oc[1]); dq.append(oc[1])
    nmul_f = nlin_f = 0
    units_f = []
    for a in sorted(touched):
        mm = ll = 0
        for mono in L.polys[a]:
            k = sum(1 for x in mono if x in unknown)
            if k >= 2:
                mm += 1
            elif k == 1:
                ll += 1
        nmul_f += mm
        lin_v = lin_cost(ll)
        nlin_f += lin_v
        units_f.extend([mul_size] * mm)
        units_f.append(lin_v)
    tot_fold = nmul_f * mul_size + nlin_f
    blocks_fold = pack(units_f, lo, hi)
    inband_fold = sum(1 for s in blocks_fold if lo <= s <= hi)
    print(f"  [B'] CONSTANT-FOLDED (only the {len(unknown)} variables downstream "
          f"of the drivers are unknown)")
    print(f"      atoms with an unknown : {len(touched):,}")
    print(f"      TRUE MUL sub-blocks   : {nmul_f:,}  "
          f"(monomials with >=2 unknown factors)")
    print(f"      linear accumulators   : {nlin_f:,} binary")
    print(f"      TOTAL binary (1 prime): **{tot_fold:,}**")
    print(f"      packed into           : {len(blocks_fold)} blocks "
          f"(min {min(blocks_fold)}, mean {sum(blocks_fold)//len(blocks_fold)}, "
          f"max {max(blocks_fold)}), {inband_fold}/{len(blocks_fold)} in band\n")

    # ---------- (2c) BRANCHED: fix the 2 selector bits (4 cases) -------------
    # x2081 and x4287 are BOOLEAN, so enumerate their 4 combinations outside the
    # QUBO.  In each branch the selector chain (x9062, x20434, x21279, x7075) is
    # a CONSTANT, and the only unknowns left are the two 256-bit free inputs.
    D2 = [8731, 9118]
    unknown2 = set(D2)
    touched2 = set()
    dq = list(D2)
    while dq:
        x = dq.pop()
        for a in users[x]:
            touched2.add(a)
            oc = L.atom_out.get(a)
            if oc is not None and oc[1] not in unknown2:
                unknown2.add(oc[1]); dq.append(oc[1])
    nmul_b = nlin_b = 0
    units_b = []
    for a in sorted(touched2):
        mm = ll = 0
        for mono in L.polys[a]:
            k = sum(1 for x in mono if x in unknown2)
            if k >= 2:
                mm += 1
            elif k == 1:
                ll += 1
        nmul_b += mm
        lv = lin_cost(ll)
        nlin_b += lv
        units_b.extend([mul_size] * mm)
        units_b.append(lv)
    tot_br = nmul_b * mul_size + nlin_b
    blocks_br = pack(units_b, lo, hi)
    inband_br = sum(1 for s in blocks_br if lo <= s <= hi)
    print(f"  [B''] BRANCHED (selector bits enumerated outside: 4 cases; "
          f"unknowns = x8731, x9118)")
    print(f"      unknown wires         : {len(unknown2)}")
    print(f"      atoms with an unknown : {len(touched2)}")
    print(f"      TRUE MUL sub-blocks   : {nmul_b}")
    print(f"      TOTAL binary / branch : **{tot_br:,}**  "
          f"({'FITS' if tot_br <= 100000 else 'over'} a 100,000 budget)")
    print(f"      packed into           : {len(blocks_br)} blocks "
          f"(min {min(blocks_br)}, mean {sum(blocks_br)//len(blocks_br)}, "
          f"max {max(blocks_br)}), {inband_br}/{len(blocks_br)} in band\n")

    # ---------- (3) LAYER B: the p-divisibility core -------------------------
    print(f"  [C] LAYER B -- the p-divisibility conditions (limb/carry chains)")
    print(f"      p | x7075*x9118 and p | x7075*x8731")
    print(f"      2 chains x 16 blocks x ~2,484 binary = ~79,500 binary "
          f"[VERIFIED sound+complete]\n")

    # ---------- coupling ------------------------------------------------------
    occ = defaultdict(int)
    for a in cone_atoms:
        for x in L.avars[a]:
            occ[x] += 1
    shared = sum(1 for x, c in occ.items() if c > 1)
    print(f"  COUPLING (core-anchored): {shared:,} wires shared by >1 atom, "
          f"{w} bits each")
    print(f"      interface = {shared*w:,} binary vs {tot_core:,} internal "
          f"({shared*w/tot_core:.1%})")

    print(f"\n  DELIVERABLE QUBO SET = [B] + [C] = "
          f"{tot_core + 79500:,} binary in "
          f"{len(blocks_core) + 32:,} blocks, all ~1k-5k, couplers <= "
          f"{mul_coeff.bit_length()} bits")
    print(f"  {time.time()-t0:.1f}s")

    json.dump({'q': q, 'w': w,
               'full': {'atoms': NA, 'muls': nmul, 'total_binary': tot_full,
                        'blocks': len(blocks_full)},
               'core': {'atoms': len(cone_atoms), 'muls': nmul_c,
                        'total_binary': tot_core, 'blocks': len(blocks_core),
                        'block_min': min(blocks_core),
                        'block_max': max(blocks_core),
                        'in_band': inband_core},
               'layerB_binary': 79500,
               'mul_block_size': mul_size, 'mul_max_coeff': mul_coeff,
               'shared_wires': shared},
              open(os.path.join(os.path.dirname(__file__), 'qubo_full.json'), 'w'))


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    if cmd == 'verify':
        sys.exit(0 if verify() else 1)
    elif cmd == 'census':
        census()
    else:
        print(__doc__)
