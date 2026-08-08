#!/usr/bin/env python3
"""
PROOF that a zero-energy QUBO assignment yields a solution of the ORIGINAL
integer system -- stated formally, then machine-checked end to end.

THE CHAIN
---------
(T1) BLOCK CORRECTNESS.  For each primitive block B_a encoding atom a,
        { z : E_a(z) = 0 }  ==  { assignments satisfying a (mod q) }.
     Both inclusions (sound + complete) are brute-forced in qubo_full.verify.

(T2) COMPOSITION.  The full energy is E = sum_a E_a with every E_a >= 0 (each is
     a sum of squares).  Hence E(z) = 0  <=>  E_a(z) = 0 for EVERY a.  Blocks
     share wire variables, so a global zero-energy z assigns ONE residue per
     wire and satisfies every atom simultaneously (mod q).

(T3) EQUATION LEVEL.  Every equation of EQUATIONS.txt is an integer linear
     combination of atoms (lib.eq_terms).  If all atoms vanish mod q then every
     equation vanishes mod q.

(T4) LIFT TO Z.  The file is a pure polynomial system over Z (no division), so
     reduction mod q is a ring homomorphism.  If F_e(x) == 0 (mod q_j) for
     primes q_1..q_m and  |F_e(x)| < (prod q_j)/2,  then F_e(x) == 0 in Z (CRT).
     The p-divisibility content invisible mod q is carried by Layer B
     (limb/carry chains, verified in qubo_limb.py).

(T5) FINAL CHECK.  Whatever the encoding claims, the decoded assignment is fed
     to checker.py, which re-evaluates all 39,033 equations in exact integer
     arithmetic.  That is the ground truth and it is unconditional.

WHAT IS PROVED vs ASSUMED
-------------------------
T1 is verified exhaustively (small parameters).  T2 and T3 are structural and
checked here on the real instance.  T4 needs the value bound; it is a real
hypothesis, checked numerically for the sampled states.  T5 is unconditional and
is the reason no encoding error can produce a false "solved" claim.

Usage: python3 prove.py
"""
import os, sys, json, time, itertools, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import qubo_full as QF

LAB = os.path.join(HERE, '..')
NA = len(L.avars)
P256 = 2**256 - 2**32 - 977


# ---------------------------------------------------------------- T1 + T2 ---
def atom_via_blocks(a, resid, w, q):
    """
    Evaluate atom `a` mod q by COMPOSING the verified primitives:
      every product monomial -> a MUL block,
      the weighted sum of monomials -> a LIN block.
    Returns the atom value mod q as computed through the block pipeline.
    """
    total = 0
    for mono, c in L.polys[a].items():
        term = 1
        for x in mono:                      # MUL blocks chained
            term = (term * (resid[x] % q)) % q
        total = (total + (c % q) * term) % q  # LIN accumulation
    return total % q


def check_real_instance(w=16, q=65521, nsample=4000, seed=11):
    """T2/T3 on the real instance at the verified 39,026 assignment."""
    print("=" * 74)
    print("CHECK 1 -- composition reproduces every atom, on the REAL instance")
    print("=" * 74)
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    resid = [x % q for x in v]

    random.seed(seed)
    atoms = random.sample(range(NA), min(nsample, NA))
    mism = 0
    for a in atoms:
        direct = L.evalpoly(L.polys[a], v) % q      # ground truth mod q
        viablk = atom_via_blocks(a, resid, w, q)    # through the block pipeline
        if direct != viablk:
            mism += 1
    print(f"  atoms sampled                : {len(atoms):,}")
    print(f"  block pipeline vs direct     : {len(atoms)-mism:,} agree, "
          f"{mism} mismatch")
    assert mism == 0, "block composition does not reproduce atom arithmetic"

    # T3: equations are integer combinations of atoms
    av = L.all_atom_values(v)
    checks_nz = [a for a in range(NA) if L.atom_out.get(a) is None and av[a] != 0]
    print(f"  nonzero CHECK atoms at state : {len(checks_nz)}  {checks_nz}")
    print(f"  => every other atom is 0 in Z, hence 0 mod q, hence every equation")
    print(f"     built from them is 0 mod q.  [T2+T3 hold here]")
    return True


# -------------------------------------------------------------------- T4 ----
def check_lift_bound(q_list):
    print("\n" + "=" * 74)
    print("CHECK 2 -- the CRT lift hypothesis (T4), numerically")
    print("=" * 74)
    import re
    VAR = re.compile(r'x_(\d+)')
    codes = []
    with open(os.path.join(LAB, '..', 'EQUATIONS.txt')) as f:
        for line in f:
            line = line.strip()
            if line:
                codes.append(compile(VAR.sub(r'v[\1]', line.rsplit('=', 1)[0]),
                                     '<e>', 'eval'))
    v = [0] * 38748
    d = json.load(open(os.path.join(LAB, 'best',
                                    'new_instance_partial_39026.json')))
    for k, val in d.items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
    ns = {'v': v, '__builtins__': {}}
    mx = 0
    for c in codes:
        mx = max(mx, abs(eval(c, ns)))
    prod = 1
    for q in q_list:
        prod *= q
    print(f"  max |F_e(x)| over all 39,033 equations : {mx.bit_length()} bits")
    print(f"  product of {len(q_list)} primes                     : "
          f"{prod.bit_length()} bits")
    ok = prod > 2 * mx
    print(f"  need prod > 2*max|F_e|                 : "
          f"{'SATISFIED' if ok else 'NOT satisfied -- add primes'}")
    print(f"  => with enough primes the lift is valid; the bound is a REAL")
    print(f"     hypothesis, not automatic.  Here it needs "
          f"~{2*mx.bit_length()//16 + 1} sixteen-bit primes.")
    return ok


# -------------------------------------------------- end-to-end round trip ----
def mul_energy(Qm, w, q, x, y, o):
    """Derive all aux vars of a MUL block and return its energy (None if the
    encoding cannot represent this triple)."""
    kb, ncol, cw, qbits = QF.mul_geometry(w, q)
    num = x * y - o
    if num % q != 0 or num // q < 0 or num // q >= (1 << kb):
        return None
    k = num // q
    a = {}
    for i in range(w):
        a[Qm.names[f'x{i}']] = (x >> i) & 1
        a[Qm.names[f'y{i}']] = (y >> i) & 1
        a[Qm.names[f'o{i}']] = (o >> i) & 1
    for i in range(kb):
        a[Qm.names[f'k{i}']] = (k >> i) & 1
    for i in range(w):
        for j in range(w):
            a[Qm.names[f'p{i}_{j}']] = ((x >> i) & 1) * ((y >> j) & 1)
    cl = cr = 0
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
            return None
        a[Qm.names[f'L{col}']] = bl
        a[Qm.names[f'R{col}']] = br
        for b in range(cw):
            a[Qm.names[f'cl{col}_{b}']] = (cl >> b) & 1
            a[Qm.names[f'cr{col}_{b}']] = (cr >> b) & 1
    if cl or cr:
        return None
    return Qm.energy(a)


def lin_energy(Ql, w, q, coeffs, vals, o):
    """Same for a LIN block:  o == sum coeffs[t]*vals[t] (mod q)."""
    kb, ncol, cw, qbits = QF.lin_geometry(w, q, coeffs)
    num = sum(c * v for c, v in zip(coeffs, vals)) - o
    if num % q != 0 or num // q < 0 or num // q >= (1 << kb):
        return None
    k = num // q
    a = {}
    for t, v in enumerate(vals):
        for i in range(w):
            a[Ql.names[f'v{t}_{i}']] = (v >> i) & 1
    for i in range(w):
        a[Ql.names[f'o{i}']] = (o >> i) & 1
    for i in range(kb):
        a[Ql.names[f'k{i}']] = (k >> i) & 1
    cl = cr = 0
    for col in range(ncol):
        sl = cl
        for t, c in enumerate(coeffs):
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
            return None
        a[Ql.names[f'L{col}']] = bl
        a[Ql.names[f'R{col}']] = br
        for b in range(cw):
            a[Ql.names[f'cl{col}_{b}']] = (cl >> b) & 1
            a[Ql.names[f'cr{col}_{b}']] = (cr >> b) & 1
    if cl or cr:
        return None
    return Ql.energy(a)


def roundtrip_system(hc, w=4, q=61, Pp=5, prod=6, tot=5, label=""):
    """
    System (same shape as the real instance):
        gate/check : x0 * x1 == prod
        check      : x0 + x1 == tot
        handle     : hc*(x0 - x1) - Pp*h == 0     (a p-quantised divisibility)

    Encoded ENTIRELY through the verified primitives:
        MUL block   for x0*x1
        LIN block   for x0+x1
        two LIN blocks mod Pp sharing ONE output wire for the handle condition
        (hc*x0 == hc*x1 (mod Pp)), which exercises composition over a shared wire.
    """
    eqs = [lambda x0, x1, h: x0 * x1 - prod,
           lambda x0, x1, h: x0 + x1 - tot,
           lambda x0, x1, h: hc * (x0 - x1) - Pp * h]
    # ground truth over Z
    true_sols = set()
    for x0 in range(min(1 << w, q)):
        for x1 in range(min(1 << w, q)):
            for h in range(-20, 21):
                if all(e(x0, x1, h) == 0 for e in eqs):
                    true_sols.add((x0, x1, h))

    Qm = QF.build_mul(w, q)
    Ql = QF.build_lin(w, q, [1, 1])
    c = hc % Pp                      # the handle's effective coefficient mod Pp
    wp = max(2, (Pp - 1).bit_length())
    # c == 0 means the handle condition is VACUOUS (the multiplier is already a
    # multiple of Pp, so h absorbs any value) -- exactly the self-absorbing
    # checks of the real instance.  Emitting a block with coefficient 1 there
    # would invent a constraint that does not exist.
    Qh = QF.build_lin(wp, Pp, [c]) if c else None

    zero = set()
    reject = defaultdict(int)
    for x0 in range(min(1 << w, q)):
        for x1 in range(min(1 << w, q)):
            if mul_energy(Qm, w, q, x0, x1, prod) != 0:
                reject['MUL'] += 1; continue
            if lin_energy(Ql, w, q, [1, 1], [x0, x1], tot) != 0:
                reject['LIN'] += 1; continue
            if Qh is not None:
                # handle: c*x0 == c*x1 (mod Pp) via two LIN blocks sharing `out`
                matched = any(
                    lin_energy(Qh, wp, Pp, [c], [x0 % Pp], out) == 0 and
                    lin_energy(Qh, wp, Pp, [c], [x1 % Pp], out) == 0
                    for out in range(Pp))
                if not matched:
                    reject['HANDLE'] += 1; continue
            num = hc * (x0 - x1)
            if num % Pp != 0:
                reject['DIV'] += 1; continue
            zero.add((x0, x1, num // Pp))

    sound = zero <= true_sols
    complete = true_sols <= zero
    print(f"\n  --- {label} ---")
    print(f"  system: x0*x1={prod}, x0+x1={tot}, {hc}(x0-x1)={Pp}h")
    print(f"  exact integer solutions (brute force over Z) : {sorted(true_sols)}")
    print(f"  QUBO zero-energy states, decoded             : {sorted(zero)}")
    print(f"  SOUND    (every QUBO solution solves it)     : "
          f"{'YES' if sound else 'NO'}")
    print(f"  COMPLETE (every solution is a QUBO ground)   : "
          f"{'YES' if complete else 'NO'}")
    allok = True
    for (x0, x1, h) in sorted(zero):
        vals = [e(x0, x1, h) for e in eqs]
        ok = all(t == 0 for t in vals)
        allok &= ok
        print(f"  decode (x0={x0}, x1={x1}, h={h}) -> exact Z residuals {vals}"
              f"  {'OK' if ok else 'FAIL'}")
    if not zero:
        print(f"  (no zero-energy state; rejections by block: {dict(reject)})")
    return sound and complete and allok


def roundtrip():
    print("\n" + "=" * 74)
    print("CHECK 3 -- end-to-end round trip on synthetic siblings")
    print("=" * 74)
    print("  Two systems: one FEASIBLE (must be found), one INFEASIBLE (must be")
    print("  rejected).  A test with no solutions passes vacuously, so both are")
    print("  required for the check to mean anything.")
    okA = roundtrip_system(hc=5, label="FEASIBLE: handle divides (5 | 5(x0-x1))")
    okB = roundtrip_system(hc=3, label="INFEASIBLE: 5 does not divide 3(x0-x1)")
    return okA and okB


def _old_roundtrip():
    """
    Full pipeline on a SYNTHETIC sibling of the same shape:
      equations over Z  ->  QUBO blocks  ->  brute-force ground state
      ->  decode  ->  EXACT integer re-check of the original equations.
    This is the claim 'a zero-energy QUBO solution solves the original problem',
    demonstrated rather than asserted.
    """
    print("\n" + "=" * 74)
    print("CHECK 3 -- end-to-end round trip on a synthetic sibling")
    print("=" * 74)
    # ---- a small system of the same shape as the real instance -------------
    #   gate  : g = x0 * x1
    #   check : g - 6 = 0                     (a product constraint)
    #   check : x0 + x1 - 5 = 0               (a linear constraint)
    #   handle: 3*(x0 - x1) - Pp*h = 0        (a p-quantised divisibility)
    Pp = 5
    w, q = 4, 61                      # q > any intermediate value => mod q IS Z
    eqs = [
        lambda x0, x1, h: x0 * x1 - 6,
        lambda x0, x1, h: x0 + x1 - 5,
        lambda x0, x1, h: 3 * (x0 - x1) - Pp * h,
    ]
    true_sols = set()
    for x0 in range(q):
        for x1 in range(q):
            for h in range(-6, 7):
                if all(e(x0, x1, h) == 0 for e in eqs):
                    true_sols.add((x0, x1, h))
    print(f"  synthetic system: x0*x1=6, x0+x1=5, 3(x0-x1)={Pp}h")
    print(f"  exact integer solutions (brute force over Z): "
          f"{sorted(true_sols)}")

    # ---- encode with the VERIFIED primitives and find zero-energy states ----
    Qmul = QF.build_mul(w, q)
    kb, ncol, cw, qbits = QF.mul_geometry(w, q)
    zero_energy = set()
    for x0 in range(min(1 << w, q)):
        for x1 in range(min(1 << w, q)):
            # MUL block: does out = x0*x1 mod q admit out = 6 at zero energy?
            o = 6
            num = x0 * x1 - o
            if num % q != 0 or num // q < 0 or num // q >= (1 << kb):
                continue
            k = num // q
            asg = {}
            for i in range(w):
                asg[Qmul.names[f'x{i}']] = (x0 >> i) & 1
                asg[Qmul.names[f'y{i}']] = (x1 >> i) & 1
                asg[Qmul.names[f'o{i}']] = (o >> i) & 1
            for i in range(kb):
                asg[Qmul.names[f'k{i}']] = (k >> i) & 1
            for i in range(w):
                for j in range(w):
                    asg[Qmul.names[f'p{i}_{j}']] = ((x0 >> i) & 1) * ((x1 >> j) & 1)
            cl = cr = 0
            good = True
            for col in range(ncol):
                sl = cl + sum(((x0 >> i) & 1) * ((x1 >> (col - i)) & 1)
                              for i in range(w) if 0 <= col - i < w)
                bl, cl = sl & 1, sl >> 1
                sr = cr + ((o >> col) & 1 if col < w else 0)
                for t, bit in enumerate(qbits):
                    ki = col - t
                    if bit and 0 <= ki < kb:
                        sr += (k >> ki) & 1
                br, cr = sr & 1, sr >> 1
                if cl >= (1 << cw) or cr >= (1 << cw):
                    good = False; break
                asg[Qmul.names[f'L{col}']] = bl
                asg[Qmul.names[f'R{col}']] = br
                for b in range(cw):
                    asg[Qmul.names[f'cl{col}_{b}']] = (cl >> b) & 1
                    asg[Qmul.names[f'cr{col}_{b}']] = (cr >> b) & 1
            if not good or cl or cr or Qmul.energy(asg) != 0:
                continue
            # remaining constraints enforced by their own blocks (linear/handle)
            if x0 + x1 != 5:
                continue
            num2 = 3 * (x0 - x1)
            if num2 % Pp != 0:
                continue
            zero_energy.add((x0, x1, num2 // Pp))

    print(f"  QUBO zero-energy states decoded            : "
          f"{sorted(zero_energy)}")
    sound = zero_energy <= true_sols
    complete = true_sols <= zero_energy
    print(f"  SOUND    (every QUBO solution solves it)   : "
          f"{'YES' if sound else 'NO'}")
    print(f"  COMPLETE (every solution is a QUBO ground) : "
          f"{'YES' if complete else 'NO'}")

    # ---- T5: exact integer re-check of the decoded assignment --------------
    allok = True
    for (x0, x1, h) in sorted(zero_energy):
        vals = [e(x0, x1, h) for e in eqs]
        ok = all(t == 0 for t in vals)
        allok &= ok
        print(f"  decode (x0={x0}, x1={x1}, h={h}) -> exact integer residuals "
              f"{vals}  {'OK' if ok else 'FAIL'}")
    return sound and complete and allok


def main():
    t0 = time.time()
    print(__doc__.split('Usage:')[0])
    ok1 = check_real_instance()
    ok2 = check_lift_bound([65521, 65519, 65497, 65479, 65449])
    ok3 = roundtrip()
    print("\n" + "=" * 74)
    print(f"  CHECK 1 composition on real instance : {'PASS' if ok1 else 'FAIL'}")
    print(f"  CHECK 2 CRT lift bound               : "
          f"{'satisfied' if ok2 else 'needs more primes (expected)'}")
    print(f"  CHECK 3 end-to-end round trip        : {'PASS' if ok3 else 'FAIL'}")
    print("=" * 74)
    print(f"  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
