"""S10 step 72: rational reconstruction and the residue battery, done properly.

The structural hypothesis about the constants is dead (see NOTEBOOK 32).  The remaining number-theoretic question is whether
the setter's residues are STRUCTURED: derived from small rationals, small
multiples of each other, or from constants elsewhere in the file.

Rational reconstruction: given r mod p, the extended-Euclid / continued-fraction
sweep on (p, r) yields the unique small a/b with a/b == r (mod p) whenever
|a|,|b| < sqrt(p/2).  If the setter built residues from small rationals, this
finds them.  It is the standard test and I had not run it.
"""
import os, sys, math, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
D0 = (v[7068] - v[2099]) % P
K2 = v[28730] % P
HUGE = 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
C1 = 33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319


def ratrec(r, m, bound=None):
    """Smallest a/b with a == b*r (mod m).  Returns (a, b, size)."""
    if bound is None:
        bound = math.isqrt(m // 2)
    a0, a1 = m, r % m
    b0, b1 = 0, 1
    while a1 > bound:
        q = a0 // a1
        a0, a1 = a1, a0 - q * a1
        b0, b1 = b1, b0 - q * b1
    return a1, b1, max(abs(a1), abs(b1))


print('=== rational reconstruction (a/b == r mod p) ===')
for nm, r in (('D0', D0), ('K2', K2), ('D0/K2', D0 * pow(K2, -1, P) % P),
              ('K2/D0', K2 * pow(D0, -1, P) % P),
              ('D0*K2', D0 * K2 % P), ('D0+K2', (D0 + K2) % P),
              ('D0-K2', (D0 - K2) % P), ('1/D0', pow(D0, -1, P)),
              ('1/K2', pow(K2, -1, P)),
              ('HUGE mod p', HUGE % P), ('C1 mod p', C1 % P)):
    a, b, sz = ratrec(r, P)
    print(f'  {nm:<12} a={str(a)[:34]:<36} b={str(b)[:22]:<24} '
          f'max digits={len(str(sz))}')

print('\n=== the two setter constants that must agree ===')
print(f'HUGE (pin 31670 on x_22152) = {HUGE}')
print(f'C1   (pin 3576  on x_6418 ) = {C1}')
print(f'HUGE - C1 = {HUGE - C1}')
print(f'(HUGE - C1) mod p = {(HUGE - C1) % P}   == D0: {(HUGE-C1) % P == D0}')
print(f'HUGE mod p = {HUGE % P}')
print(f'C1   mod p = {C1 % P}')
print(f'HUGE // p  = {HUGE // P}    C1 // p = {C1 // P}')
print(f'gcd(HUGE, C1) = {math.gcd(HUGE, C1)}')
for k in range(2, 200):
    if (HUGE - k * C1) % P == 0:
        print(f'  *** HUGE == {k}*C1 (mod p)')
    if (C1 - k * HUGE) % P == 0:
        print(f'  *** C1 == {k}*HUGE (mod p)')

# ---- residues vs constants > p, done properly ------------------------------
consts = set()
for a in range(L.NA):
    for m, c in L.polys[a].items():
        if abs(c) > P:
            consts.add(abs(c))
print(f'\n=== constants strictly larger than p: {len(consts)} ===')
small_res = [c for c in consts if c % P < 2**80]
print(f'  with (c mod p) < 2^80 : {len(small_res)}   '
      f'(random expectation ~{len(consts) * 2**80 / P:.2e})')
if small_res:
    for c in sorted(small_res)[:5]:
        print(f'    c={c}  ->  c mod p = {c % P}')
quo = [c // P for c in consts]
print(f'  quotient c//p : min {min(quo)} max {max(quo)} '
      f'({len(str(min(quo)))}-{len(str(max(quo)))} digits)')

# do any two constants share a residue mod p?
byres = collections.Counter(c % P for c in consts)
dup = [r for r, n in byres.items() if n > 1]
print(f'  distinct residues: {len(byres)} of {len(consts)}; repeats: {len(dup)}')

# ---- the failing equation values -------------------------------------------
print('\n=== the seven failing equation values ===')
fail = L.failing_eqs(av)


def smallfac(n, lim=200000):
    out = []
    n = abs(n)
    for q in range(2, lim):
        if q * q > n: break
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q; e += 1
            out.append((q, e))
    return out, n


g = 0
for e in fail:
    val = L.eq_value(e, av)
    g = math.gcd(g, abs(val))
    f, rem = smallfac(val)
    print(f'  eq {e:<6} {len(str(abs(val)))} digits  small factors {f[:6]}')
print(f'  gcd of the seven residual values: {g}  ({len(str(g))} digits)')
if g > 1:
    fg, remg = smallfac(g)
    print(f'    factors of the gcd: {fg[:10]}  remaining {len(str(remg))} digits')
    print(f'    is the gcd divisible by p? {g % P == 0}')
