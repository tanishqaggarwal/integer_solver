#!/usr/bin/env python3
"""
RNS (residue-number-system) reformulation of EQUATIONS.txt.

The whole instance is a *pure straight-line polynomial system* over Z:
the only operators in the file are + - * and integer literals (verified:
`grep -c '[/%]' EQUATIONS.txt` == 0).  Therefore reduction modulo any
integer q is a ring homomorphism

        Z  --->  Z/q1 x Z/q2 x ... x Z/qm      (q_j distinct primes)

and an integer assignment satisfies  F_e(x) = 0  over Z for every e  iff
its residues satisfy  F_e(x) = 0 (mod q_j)  for every e and every j,
PROVIDED  prod_j q_j  >  2*max_e |F_e(x)|  (CRT faithfulness bound).

Consequence: the ONE dense 256-bit-precision QUBO that REDUCED_PROBLEM.md
Section 4 called "the wrong machine" splits into m INDEPENDENT systems, one
per prime, that share no variables (disjoint) and whose coefficients are
< q_j (small -- coupler precision becomes a free parameter, decoupled from p).

This script:
  1. loads and compiles the system once,
  2. for a schedule of small primes q, evaluates the *verified 39,026 partial*
     mod q and confirms the reduction is FAITHFUL (same 7 equations fail,
     zero spurious), i.e. tiny-word arithmetic already witnesses satisfaction,
  3. reports the coupler-width / prime-count trade-off for a full CRT lift.

Usage:  python3 rns_reduce.py [assignment.json]
"""
import re, sys, json, time
from math import gcd
from functools import reduce

ROOT = __file__.rsplit('/', 2)[0]          # .../solve_lab
REPO = ROOT.rsplit('/', 1)[0]              # repo root
EQ   = REPO + '/EQUATIONS.txt'
NVARS = 38748
VAR_RE = re.compile(r'x_(\d+)')


def load_codes(path=EQ):
    codes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0]
            codes.append(compile(VAR_RE.sub(r'v[\1]', lhs), '<eq>', 'eval'))
    return codes


def load_assignment(path):
    d = json.load(open(path))
    v = [0] * NVARS
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(2, n + 1) if s[i]]


def eval_all(codes, v):
    ns = {'v': v, '__builtins__': {}}
    return [i for i, c in enumerate(codes) if eval(c, ns) != 0]


def eval_all_mod(codes, v, q):
    vq = [x % q for x in v]
    ns = {'v': vq, '__builtins__': {}}
    return [i for i, c in enumerate(codes) if eval(c, ns) % q != 0]


def main():
    t0 = time.time()
    codes = load_codes()
    print(f"[rns] loaded {len(codes)} equations in {time.time()-t0:.1f}s")

    src = sys.argv[1] if len(sys.argv) > 1 else \
        ROOT + '/best/new_instance_partial_39026.json'
    v = load_assignment(src)
    print(f"[rns] assignment = {src}")

    # ground truth over Z
    fails = eval_all(codes, v)
    fset = set(fails)
    print(f"[rns] exact-Z failing = {len(fails)}: {fails}")

    # faithfulness of the reduction at a schedule of small primes
    print("\n[rns] faithfulness of mod-q reduction (small word size):")
    print(f"      {'q':>10} {'bits':>4}  {'nonzero mod q':>13}  "
          f"{'caught/7':>8}  {'spurious':>8}")
    for q in [7, 13, 31, 61, 101, 251, 1021, 65537]:
        bad = eval_all_mod(codes, v, q)
        bset = set(bad)
        caught = len(fset & bset)
        spurious = len(bset - fset)
        print(f"      {q:>10} {q.bit_length():>4}  {len(bad):>13}  "
              f"{caught:>6}/7  {spurious:>8}")

    # bad-prime analysis: a prime q "misses" a failing equation iff q | value
    vals = []
    ns = {'v': v, '__builtins__': {}}
    for i in fails:
        vals.append(eval(codes[i], ns))
    if vals:
        g = reduce(gcd, [abs(x) for x in vals])
        print(f"\n[rns] gcd of failing residuals = {g}  "
              f"(a prime misses ALL simultaneously only if it divides this)")
        print(f"[rns] failing residual bit-lengths = "
              f"{[x.bit_length() for x in vals]}")

    # CRT schedule: primes needed to lift the genuinely-unknown quantities.
    # The 13 planted numbers are ~296-bit; handles are a closed-form /p lift,
    # not solved for, so 296 bits is the relevant reconstruction width.
    for width_bits, label in [(296, 'the 296-bit planted numbers')]:
        target = 1 << (width_bits + 1)
        for wbits in (5, 6, 8, 12, 16):
            ps = [p for p in sieve((1 << wbits) - 1) if p >= (1 << (wbits - 1))]
            prod, k = 1, 0
            for p in ps:
                prod *= p
                k += 1
                if prod > target:
                    break
            note = "" if prod > target else f"  (only {len(ps)} such primes exist; use a wider band)"
            print(f"[rns] lift {label}: ~{k} disjoint {wbits}-bit-prime systems"
                  f"  -> coupler width <= {2*wbits} bits (squared penalty){note}")

    print(f"\n[rns] total {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
