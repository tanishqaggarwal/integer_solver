"""S11 step 10: number theory on the BROADCAST constants (not derived residues).

Earlier sessions did structure tests on D0 and K2, which are derived.  These
four are the instance's own loaded constants -- the ones the circuit asserts
against -- and they have never been examined.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
cnt = collections.Counter(v)
C = {}
for val, n in cnt.most_common(8):
    if n >= 40 and val not in (0, 1): C[n] = val
names = {}
for n, val in C.items(): names[n] = val
print('constants:', {n: str(v)[:20]+'...' for n, v in names.items()})

def leg(a, m=P): return pow(a % m, (m - 1) // 2, m)
def is_x(x):
    y2 = (pow(x, 3, P) + 7) % P
    return leg(y2) == 1
    r = val % P
    print(f'  x{n:<4} r = {str(r)[:26]}...  QR: {is_x(r)}  '
          f'r == Gx? {r == GX}  r == n? {r == N}')
print('\n=== pairwise relations mod p ===')
ks = sorted(names)
for i in range(len(ks)):
    for j in range(len(ks)):
        if i == j: continue
        a, b = names[ks[i]] % P, names[ks[j]] % P
        if b == 0: continue
        q = a * pow(b, -1, P) % P
        small = q if q < 10**6 else (q - P if P - q < 10**6 else None)
        if small is not None:
            print(f'  x{ks[i]} / x{ks[j]} = {small} (mod p)   *** SMALL')
        # continued-fraction / rational reconstruction of the ratio
        num, den, r0, r1 = 1, 0, P, q
        s0, s1 = 0, 1
        lim = 1 << 128
        while r1 > lim:
            qq = r0 // r1
            r0, r1 = r1, r0 - qq * r1
            s0, s1 = s1, s0 - qq * s1
        if abs(s1) < 10**8 and r1 < 10**8:
            print(f'  x{ks[i]} / x{ks[j]} = {r1}/{s1} (mod p)  *** RATIONAL')

print('\n=== additive / power relations ===')
for n, val in names.items():
    r = val % P
    for lbl, t in [('r^2', pow(r, 2, P)), ('r^3', pow(r, 3, P)),
                   ('1/r', pow(r, -1, P) if r else 0), ('r^3+7', (pow(r,3,P)+7)%P)]:
        for n2, val2 in names.items():
            if t == val2 % P: print(f'  x{n} {lbl} == x{n2}  ***')
    for k in range(2, 64):
        if (r * k) % P in {u % P for u in names.values()}:
            print(f'  {k} * x{n} matches another constant  ***')
print('\n=== 296-bit constants: byte structure ===')
for n, val in names.items():
    if val.bit_length() == 296:
        b = val.to_bytes(37, 'big')
        print(f'  x{n}: {b.hex()}')
        print(f'      first 5 bytes {b[:5].hex()}  last 32 {b[5:].hex()}')
        top, bot = val >> 256, val & ((1 << 256) - 1)
        print(f'      top40 {top} ({hex(top)})  bottom256 QR? {is_x(bot % P)}')
