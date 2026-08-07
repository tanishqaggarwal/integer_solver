#!/usr/bin/env python3
"""Agent Z: independent audit of agent AB's THEOREM D (generic-group bound for the
weight predicate).  Three things to check, as tasked:
  (1) is min(|D0|,|D1|) the right normaliser?
  (2) does the single-root argument survive the automorphism group {+-1, +-lambda, +-lambda^2}?
  (3) does  B = 128  =>  m >= 2^127.5  actually follow?
Class sizes are recomputed by an exact digit-DP over k in [0,N), NOT over [0,2^256).
"""
from math import comb, log2, sqrt

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
def L(x): return log2(x) if x > 0 else float('-inf')

# ---------- exact count of k in [0,N) with popcount(k) <= B  (digit DP) ----------
def count_le(B):
    if B < 0:
        return 0
    tot = 0
    ones = 0
    for i in range(255, -1, -1):
        if (N >> i) & 1:
            r = B - ones
            if r >= 0:
                tot += sum(comb(i, t) for t in range(0, min(i, r) + 1))
            ones += 1
            if ones > B:
                break
    return tot

assert count_le(256) == N, "digit-DP does not reproduce N"
print("digit-DP self-check: #{k<N} = N  ->", count_le(256) == N)
print("max popcount over k<N :", max(B for B in range(257) if count_le(B) > count_le(B - 1)),
      "  (matches AB's 'free unconditional w <= 255')")

print()
print("=" * 92)
print("(1) IS min(|D0|,|D1|) THE RIGHT NORMALISER?   -- yes, with a factor-2 correction")
print("=" * 92)
print("""  Shoup simulation: run the game with k a formal variable.  The simulated view is
  independent of k.  Real and simulated views diverge only if some pair of held elements
  sigma(a_i + b_i k), sigma(a_j + b_j k) collides, i.e. (b_i - b_j) k = a_j - a_i.  N is
  PRIME, so each pair with b_i != b_j has EXACTLY ONE root; pairs with b_i = b_j, a_i != a_j
  have none.  So Bad = {k : the view differs from the simulation} has |Bad| <= C(m,2) and --
  crucially -- Bad is fixed by the simulation, independent of k.  Then for D_b uniform on a
  set of size |D_b|:
        |Pr_{D_b}[A=1] - Pr[A_sim=1]|  <=  |Bad| / |D_b|
  and therefore
        Adv = |Pr_{D0}[A=1] - Pr_{D1}[A=1]|  <=  |Bad|*(1/|D0| + 1/|D1|)
                                             <=  (m^2/2)*(1/|D0| + 1/|D1|)
                                             <=  m^2 / min(|D0|,|D1|).
  So the NORMALISER min(|D0|,|D1|) is CORRECT -- it is the dominant term of the harmonic sum.
  AB's CONSTANT is not: AB writes Adv <= m^2 / (2 min), which drops the second side entirely
  and is a factor 2 TIGHTER than the argument supports.  Direction of the error: it OVERSTATES
  the barrier, by exactly 0.5 bits in m.""")

print()
print("=" * 92)
print("(2) DOES THE SINGLE-ROOT ARGUMENT SURVIVE THE AUTOMORPHISM GROUP?")
print("=" * 92)
print("""  Two separate questions, and they have different answers.

  (a) Can lambda or negation break the AFFINE FORM?  NO.  Both are multiplication by a fixed
      scalar of Z_N (lambda has order 3 in Z_N*, negation is -1).  Applying either to
      sigma(a + b k) yields sigma(lambda*a + lambda*b*k) -- still of the form sigma(a'+b'k).
      So every element an algorithm can hold is still affine in k, and every pairwise
      collision is still ONE equation over the field Z_N with exactly ONE root.
      => the single-root argument SURVIVES.  This is the part AB gets right.

  (b) Can the ENCODING create extra collisions?  YES, and AB's knob line ('coordinate encoding
      excluded') is doing real work here.  If the encoding is by x-coordinate, sigma(x) and
      sigma(-x) are the SAME string; with GLV the whole order-6 orbit {+-1,+-lam,+-lam^2}
      collapses.  Then a pair collides if a_i + b_i k = u*(a_j + b_j k) for ANY u in a group
      of size AUT, i.e. AUT affine equations per pair -- each still with one root, but
      |Bad| <= AUT * C(m,2).  The bound degrades by sqrt(AUT).""")
for aut, nm in ((1, "opaque encoding (AB's model)"), (2, "x-coordinate only (negation)"),
                (6, "x-coordinate + GLV (the real secp256k1 attack surface)")):
    print("     AUT = %d  (%-42s) : m lower bound degrades by 2^%.2f" % (aut, nm, 0.5 * log2(aut)))

print()
print("=" * 92)
print("(3) DOES  B = 128  =>  m >= 2^127.5  FOLLOW?")
print("=" * 92)
print("   B    |{w<=B}| (k<N)   |{w>B}| (k<N)     min        AB: sqrt(min)   corrected sqrt(min/2)"
      "   +AUT=6")
for B in (10, 20, 24, 30, 40, 56, 100, 128, 148, 152, 198, 245, 255):
    a = count_le(B); b = N - a
    mn = min(a, b)
    ab = 0.5 * L(mn)                       # AB: m >= sqrt(2*eps*min), eps=1/2
    zz = 0.5 * L(mn) - 0.5                 # corrected constant: m >= sqrt(eps*min), eps=1/2
    au = zz - 0.5 * log2(6)                # and with the order-6 automorphism group
    print("  %4d   2^%7.2f       2^%7.2f      2^%7.2f     m>=2^%6.2f    m>=2^%6.2f          m>=2^%6.2f"
          % (B, L(a), L(b), L(mn), ab, zz, au))

print()
print("  AB's arithmetic reproduces: for B=128, min = 2^%.2f and sqrt(min) = 2^%.2f -> AB's 2^127.5"
      % (L(min(count_le(128), N - count_le(128))), 0.5 * L(min(count_le(128), N - count_le(128)))))
print("  So (3) FOLLOWS from AB's stated inequality.  But the inequality's constant is a factor 2")
print("  too tight (item 1), so the defensible figure is m >= 2^127.0, and once the automorphism")
print("  group is admitted into the model (item 2b) it is m >= 2^125.7.")

print()
print("=" * 92)
print("THE ONE SUBSTANTIVE PROBLEM WITH AB's HEADLINE")
print("=" * 92)
print("""  AB writes: 'B=128 needs m >= 2^127.5 (solving costs 2^126.5) -- DECIDING THE WEIGHT
  PREDICATE IS AS HARD AS SOLVING'.  Taken literally the numbers say deciding is HARDER than
  solving, which cannot be true: any algorithm that solves also decides.  The contradiction is
  not in the mathematics, it is a MODEL MISMATCH -- 2^127.5 is a generic-group bound computed
  with the automorphisms EXCLUDED, while 2^126.5 is a concrete rho cost computed with the
  automorphisms INCLUDED (the sqrt(6) of negation + GLV).  Put both in the same model:""")
mn128 = min(count_le(128), N - count_le(128))
print("     GGM lower bound, AB's constant, no automorphisms : m >= 2^%.2f" % (0.5 * L(mn128)))
print("     GGM lower bound, corrected constant, AUT = 6     : m >= 2^%.2f"
      % (0.5 * L(mn128) - 0.5 - 0.5 * log2(6)))
print("     concrete rho with negation + GLV                 : 2^126.5")
print("  -> 2^125.7 <= 2^126.5.  The inequality points the RIGHT way once the models agree, and")
print("     AB's qualitative conclusion is unharmed: deciding w <= 128 admits NO generic")
print("     shortcut and costs the same as solving to within the same sqrt(6) factor.")
print("     What does not survive is the strict claim that deciding is harder than solving.")

print()
print("=" * 92)
print("CROSS-CHECK: does the corrected MITM meet the generic bound for SMALL classes?")
print("=" * 92)
V128 = [0] * 129
s = 0
for j in range(129):
    s += comb(128, j); V128[j] = s
print("   B    GGM LB (corr.)   corrected MITM   gap (bits)")
for B in (10, 20, 24, 30, 40, 56):
    mn = min(count_le(B), N - count_le(B))
    lb = 0.5 * L(mn) - 0.5
    mitm = L(sqrt(max(B, 1)) * V128[min((B + 1) // 2, 128)])
    print("  %4d     2^%6.2f        2^%6.2f       %+5.2f" % (B, lb, mitm, mitm - lb))
print("  -> the corrected MITM is within ~1.5 bits of the generic lower bound across the whole")
print("     actionable range.  AB's claim 'within 2^1 of optimal at B=20' is CONFIRMED.")
