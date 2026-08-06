"""S10 step 71: the number theory I should have done first.

p = 2^256 - 2^32 - 977 is the secp256k1 FIELD prime.  If this instance is a curve
computation then the binding residues are not random 256-bit numbers -- they are
coordinates, and the trapdoor is point arithmetic rather than a covering design.

Tests, all cheap and all decisive:
  * is (D0, K2) a point on y^2 = x^3 + 7 ?
  * is D0 (or K2, or any constant mod p) a valid x-coordinate?
  * do the group order n, generator G, or their multiples appear as constants?
  * are the 7-digit multipliers prime / structured?
"""
import os, sys, collections, math, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141   # secp256k1 order
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
B = 7

print(f'p (secp256k1 field prime) = {P}')
print(f'n (group order)           = {N}')
print(f'p == 2^256-2^32-977 : {P == 2**256 - 2**32 - 977}')
print(f'p mod 4 = {P % 4}  (3 => sqrt is a^((p+1)/4))')

v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
D0 = (v[7068] - v[2099]) % P
K2 = v[28730] % P
print(f'\nD0 = {D0}')
print(f'K2 = {K2}')


def on_curve(x, y):
    return (y * y - (x * x * x + B)) % P == 0


def is_x(x):
    r = (x * x * x + B) % P
    return pow(r, (P - 1) // 2, P) == 1


def sqrt_p(a):
    r = pow(a, (P + 1) // 4, P)
    return r if r * r % P == a % P else None


print('\n=== is (D0, K2) a point on y^2 = x^3 + 7 ? ===')
for nm, (x, y) in (('(D0,K2)', (D0, K2)), ('(K2,D0)', (K2, D0))):
    print(f'  {nm}: on curve = {on_curve(x, y)}')
print(f'  D0 a valid x-coordinate: {is_x(D0)}')
print(f'  K2 a valid x-coordinate: {is_x(K2)}')
if is_x(D0):
    print(f'     y(D0) = {sqrt_p((D0**3+B) % P)}')
if is_x(K2):
    print(f'     y(K2) = {sqrt_p((K2**3+B) % P)}')

print('\n=== do n, G_x, G_y appear among the instance constants? ===')
consts = set()
for a in range(L.NA):
    for m, c in L.polys[a].items():
        if abs(c) > 10**9:
            consts.add(abs(c))
print(f'large constants collected: {len(consts)}')
special = {'p': P, 'n': N, 'G_x': GX, 'G_y': GY, 'p-n': P - N,
           '(p-1)/2': (P - 1) // 2, 'p+1': P + 1, 'n-1': N - 1, 'b=7': B}
for nm, val in special.items():
    hits = [c for c in consts if c == val]
    print(f'  {nm:<8} present as a literal: {bool(hits)}')
    modhits = [c for c in consts if c % P == val % P and c != val]
    if modhits:
        print(f'      but {len(modhits)} constants are == {nm} (mod p): '
              f'{sorted(modhits)[:2]}')

print('\n=== constants mod p as x-coordinates ===')
xok = 0
tested = 0
for c in sorted(consts):
    r = c % P
    if r == 0:
        continue
    tested += 1
    if is_x(r):
        xok += 1
print(f'  {xok} of {tested} constants have (c mod p)^3 + 7 a QR '
      f'(expected ~{tested//2} if random)')

print('\n=== also: n as a modulus? do constants reduce nicely mod n? ===')
small_mod_n = [c for c in consts if c % N < 2**64]
print(f'  constants with (c mod n) < 2^64: {len(small_mod_n)}')
small_mod_p = [c for c in consts if c % P < 2**64]
print(f'  constants with (c mod p) < 2^64: {len(small_mod_p)}')

print('\n=== the 7-digit multipliers ===')
mults = sorted(set(abs(c) for a in range(L.NA) for m, c in L.polys[a].items()
                   if 10**6 < abs(c) < 10**8))
print(f'  count: {len(mults)}   sample: {mults[:12]}')


def isprime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True


pr = [m for m in mults if isprime(m)]
print(f'  prime among them: {len(pr)} of {len(mults)} '
      f'(expected ~{int(len(mults)/17)} if random 7-digit)')
g = 0
for m in mults: g = math.gcd(g, m)
print(f'  gcd of all multipliers: {g}')
