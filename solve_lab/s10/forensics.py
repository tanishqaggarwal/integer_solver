"""S10 step 67: NUMBER-THEORETIC FORENSICS on the setter's constants.

Every attack so far treated the HUGE constants as opaque. They are not opaque:
they are the ONLY source of values in the instance that are not multiples of p,
so they are the trapdoor's seed.  If the setter generated them with structure,
the arithmetic will show it.

Collect every large integer literal in the instance and test:
  * size vs p, quotient and residue mod p
  * gcds, pairwise differences and ratios mod p
  * whether the two binding residues D0, K2 are related by a small multiplier
  * Legendre symbols
  * whether any constant is a multiple of p plus a small/structured remainder
"""
import os, sys, collections, math, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))

# ---- collect every large constant appearing in any atom ---------------------
consts = collections.Counter()
for a in range(L.NA):
    for m, c in L.polys[a].items():
        if abs(c) > 10**30:
            consts[abs(c)] += 1
print(f'distinct large (>10^30) constants in atoms: {len(consts)}')
big = sorted(consts)
for c in big[:6]:
    print(f'   {c}  (x{consts[c]})  {len(str(c))} digits')

print(f'\np = {P}   ({len(str(P))} digits, {P.bit_length()} bits)')
print(f'\n{"constant":>12} {"bits":>5} {"q = c//p":>22} {"c mod p (first 30)":>32}')
info = []
for c in big:
    q, r = divmod(c, P)
    info.append((c, q, r))
    print(f'{str(c)[:12]:>12} {c.bit_length():>5} {str(q)[:22]:>22} {str(r)[:30]:>32}')

# ---- is any constant a clean multiple of p, or p times something small? -----
print('\n=== structure of the quotients c // p ===')
for c, q, r in info:
    print(f'  c={str(c)[:16]}...  q bits={q.bit_length():<4} r bits={r.bit_length():<4} '
          f'p|c: {r == 0}')

# ---- pairwise gcds and differences -----------------------------------------
print('\n=== pairwise gcd of the constants ===')
g = 0
for c in big:
    g = math.gcd(g, c)
print(f'gcd of all large constants = {g}')
print('gcd of pairwise differences:')
gd = 0
for i in range(len(big)):
    for j in range(i + 1, len(big)):
        gd = math.gcd(gd, abs(big[i] - big[j]))
print(f'  {gd}   ({len(str(gd))} digits)' if gd else '  0')

# ---- the two binding residues ----------------------------------------------
D0 = (v[7068] - v[2099]) % P
K2 = v[28730] % P
print(f'\n=== the two binding residues ===')
print(f'D0 = {D0}')
print(f'K2 = {K2}')
print(f'D0 bits {D0.bit_length()}, K2 bits {K2.bit_length()}')
print(f'Legendre(D0) = {pow(D0,(P-1)//2,P) == 1}')
print(f'Legendre(K2) = {pow(K2,(P-1)//2,P) == 1}')
ratio = D0 * pow(K2, -1, P) % P
print(f'D0/K2 mod p = {ratio}   ({ratio.bit_length()} bits)')
inv = K2 * pow(D0, -1, P) % P
print(f'K2/D0 mod p = {inv}   ({inv.bit_length()} bits)')
print(f'D0 + K2 mod p = {(D0+K2) % P}')
print(f'D0 - K2 mod p = {(D0-K2) % P}')
print(f'D0 * K2 mod p = {D0*K2 % P}')
for k in range(2, 60):
    if (D0 - k * K2) % P == 0:
        print(f'  *** D0 == {k} * K2 (mod p)')
    if (K2 - k * D0) % P == 0:
        print(f'  *** K2 == {k} * D0 (mod p)')
    if (D0 * k - K2) % P == 0:
        print(f'  *** {k}*D0 == K2 (mod p)')

# ---- are the residues related to the constants? ----------------------------
print('\n=== residues vs the constants (mod p) ===')
for c, q, r in info:
    for name, val in (('D0', D0), ('K2', K2)):
        if r == val:
            print(f'  *** {name} == constant {str(c)[:20]}... mod p')
        d = (r - val) % P
        if d.bit_length() < 60:
            print(f'  {name} differs from c={str(c)[:16]}... mod p by a '
                  f'{d.bit_length()}-bit number: {d}')

# ---- small-factor forensics -------------------------------------------------
print('\n=== small prime factors ===')
def smallfac(n, lim=100000):
    out = []
    for q in range(2, lim):
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q; e += 1
            out.append((q, e))
        if q * q > n:
            break
    return out, n

for name, val in (('D0', D0), ('K2', K2), ('D0-K2', (D0 - K2) % P),
                  ('D0+K2', (D0 + K2) % P)):
    f, rem = smallfac(val)
    print(f'  {name:<8} small factors {f[:8]}  remaining {len(str(rem))} digits')
for c in big[:8]:
    f, rem = smallfac(c)
    print(f'  const {str(c)[:14]}.. small factors {f[:8]}  remaining {len(str(rem))} digits')
