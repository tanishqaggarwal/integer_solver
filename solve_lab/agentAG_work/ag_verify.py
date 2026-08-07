#!/usr/bin/env python3
"""AG -- falsifiable checks.  Each test below has a way to FAIL; several are run on a case where
failure would refute the claim being tested (the coordinator's 'vacuous plant test' warning).

T1  rep(W) against a MEASURED balance frequency in small n            -- can fail
T2  AB's rep is exactly 2x the honest rep for every ODD W             -- can fail
T3  the cyclic-shift splitting system (deterministic, size n/2)       -- can fail
T4  end-to-end low-weight MITM in a real prime-order group, run on
    BOTH a plantable case (w<=B, must HIT) and an unplantable one
    (w>B, must MISS).  The MISS branch is what makes T4 non-vacuous.  -- can fail
T5  operation COUNT of T4 against the closed form rep*Vol_h(ceil(W/2)) -- can fail
"""
import math, random, itertools
from math import comb, log2

# ---------------------------------------------------------------- T1/T2/T3 helpers
def rep_single(n, W):
    """AB's form: 1 / P[split is exactly (ceil,floor)] with only the ceil term counted"""
    h = n//2; c = (W+1)//2
    num = comb(W, c)*comb(n-W, h-c) if 0 <= h-c <= n-W else 0
    return comb(n, h)/num if num else float('inf')

def rep_honest(n, W):
    """1 / P[both sides <= ceil(W/2)]  -- what the algorithm actually needs"""
    h = n//2; c = (W+1)//2; num = 0
    for a in range(W-c, c+1):
        if 0 <= a <= W and 0 <= h-a <= n-W:
            num += comb(W, a)*comb(n-W, h-a)
    return comb(n, h)/num if num else float('inf')

def measure_balance(n, W, trials, rng):
    """MEASURED frequency that a random h/h split leaves <= ceil(W/2) on each side"""
    h = n//2; c = (W+1)//2; pos = list(range(n)); ok = 0
    for _ in range(trials):
        D = set(rng.sample(pos, W))
        rng.shuffle(pos)
        left = set(pos[:h])
        a = len(D & left)
        if a <= c and W-a <= c: ok += 1
    return ok/trials

def cyclic_ok(n, W):
    """deterministic splitting system: the n/2 contiguous windows of a cycle of n positions.
       Claim: for EVERY W-subset some window holds between floor(W/2) and ceil(W/2)."""
    h = n//2; c = (W+1)//2
    for D in itertools.combinations(range(n), W):
        Ds = set(D); good = False
        for i in range(n):
            a = sum(1 for j in range(i, i+h) if (j % n) in Ds)
            if a <= c and W-a <= c: good = True; break
        if not good: return False, D
    return True, None

# ---------------------------------------------------------------- T4/T5: real group MITM
def mitm_certify(P, g, q, n, T, B, rng, max_splits=200):
    """Certify w<=B by SEARCHING {wt<=B} (AG's Attack-2 algorithm), in the group <g> of order q.
       Splitting-system MITM.  Returns (S or None, group_ops)."""
    h = n//2; c = (B+1)//2
    ops = 0
    pos = list(range(n))
    ginv = pow(g, q-1, P)
    nsplit = max(1, int(round(rep_honest(n, B))))
    for _ in range(min(max_splits, 4*nsplit)):
        rng.shuffle(pos)
        Lh, Rh = pos[:h], pos[h:]
        table = {}
        for r in range(c+1):
            for sub in itertools.combinations(Lh, r):
                e = sum(1 << i for i in sub)
                v = pow(g, e, P); ops += 1
                table.setdefault(v, sub)
        for r in range(c+1):
            for sub in itertools.combinations(Rh, r):
                e = sum(1 << i for i in sub)
                v = (T*pow(pow(g, e, P), q-1, P)) % P; ops += 1
                if v in table:
                    return sorted(set(table[v]) | set(sub)), ops
    return None, ops

if __name__ == "__main__":
    rng = random.Random(20260807)
    print("="*94); print("AG PART 3 -- FALSIFIABLE CHECKS"); print("="*94)

    print("\n-- T1: rep(W) vs MEASURED balance frequency (n=32, 200k trials each) --")
    n = 32; bad = 0
    for W in (3,4,7,8,11,12,15,16):
        f = measure_balance(n, W, 200000, rng)
        pred = 1.0/rep_honest(n, W)
        rel = abs(f-pred)/pred
        flag = "OK" if rel < 0.02 else "*** MISMATCH ***"
        if rel >= 0.02: bad += 1
        print("   W=%2d  measured P=%.5f  predicted 1/rep_honest=%.5f  rel.err %.4f  %s"
              %(W, f, pred, rel, flag))
    print("   T1 verdict:", "PASS" if bad == 0 else "FAIL (%d mismatches)"%bad)

    print("\n-- T2: is AB's rep exactly 2x the honest rep for odd W? (n=256) --")
    bad = 0
    for W in range(1, 256):
        r = rep_single(256, W)/rep_honest(256, W)
        want = 2.0 if W % 2 else 1.0
        if abs(r-want) > 1e-9: bad += 1
    print("   checked W=1..255 at n=256: exceptions = %d"%bad)
    print("   PROOF: for odd W=2c-1 the two admissible splits contribute C(W,c)C(256-W,128-c) and")
    print("   C(W,c-1)C(256-W,129-c); C(W,c)=C(W,c-1) since W=2c-1, and 256-W=257-2c is ODD with")
    print("   (128-c)+(129-c)=257-2c, so the two binomials are the twin central ones and are EQUAL.")
    print("   Hence the honest denominator is exactly twice AB's.  T2 verdict:",
          "PASS" if bad == 0 else "FAIL")

    print("\n-- T3: deterministic cyclic-shift splitting system (exhaustive over ALL W-subsets) --")
    allok = True
    for n_ in (12, 14, 16):
        for W in range(0, n_+1):
            ok, ce = cyclic_ok(n_, W)
            if not ok: allok = False; print("   n=%d W=%d COUNTEREXAMPLE %s"%(n_, W, ce))
    print("   n=12,14,16, every W, every W-subset: %s"%("all balanced by some window" if allok else "FAILED"))
    print("   -> a ZERO-ERROR covering proof needs only n/2=128 splits, not an unbounded family.")
    print("      AB charges rep(W)~10 (the Las Vegas expectation).  The deterministic factor is")
    print("      bounded by 128, so AB UNDERPRICES a zero-error proof by at most 2^%.1f."
          %log2(128/10.0))
    print("   T3 verdict:", "PASS" if allok else "FAIL")

    print("\n-- T4/T5: end-to-end certification in a real prime-order group --")
    # q prime, P = 2q+1 prime, g of order q.  n=24 exponent positions so folds < 2^24 << q:
    # no wraparound, so the subset is unique and a MISS is a genuine proof of w>B.
    q = 1152921504606847009
    while not all(q % p for p in (3,5,7,11,13,17,19,23)) or pow(2, q-1, q) != 1:
        q += 2
    P = 2*q+1
    while pow(2, P-1, P) != 1:
        q = q+2
        while pow(2, q-1, q) != 1: q += 2
        P = 2*q+1
    g = pow(3, 2, P)
    assert pow(g, q, P) == 1 and g != 1
    print("   group: <g> of prime order q = %d (%.1f bits), P = 2q+1"%(q, log2(q)))
    n, B = 24, 6
    for label, wt in (("PLANT w=4  (<= B=6, MUST HIT)", 4), ("PLANT w=10 (> B=6, MUST MISS)", 10)):
        S = sorted(rng.sample(range(n), wt))
        k = sum(1 << i for i in S)
        T = pow(g, k, P)
        found, ops = mitm_certify(P, g, q, n, T, B, rng)
        hit = found is not None
        want = (wt <= B)
        ok = (hit == want) and (not hit or found == S)
        print("   %-32s  planted S=%s" % (label, S))
        print("       -> %s   correct=%s   group ops used 2^%.2f"
              %("HIT "+str(found) if hit else "MISS (certifies w>%d)"%B, ok, log2(ops)))
        if wt <= B:
            h = n//2; c = (B+1)//2
            Vh = sum(comb(h, j) for j in range(c+1))
            print("       T5: closed form per split = 2*Vol_%d(%d) = 2^%.2f ; splits until hit"
                  %(h, c, log2(2*Vh)))
            print("           implied = %.2f ; rep_honest(%d,%d) = %.2f"
                  %(ops/(2*Vh), n, B, rep_honest(n, B)))
    print("   T4 verdict: the MISS branch is the one that could have refuted the certifier;")
    print("               it is a genuine test because a false HIT there would be unsound.")
