import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
cnt = collections.Counter(v)
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141   # secp256k1 order
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
known = {'p_secp256k1': P, 'n_secp256k1': N, 'Gx': GX, 'Gy': GY,
         '2^256': 2**256, '2^255': 2**255, 'p-1': P-1, 'n-1': N-1}
print(f'{"mult":>6}  value')
for val, n in cnt.most_common(14):
    if n < 5 or val in (0, 1): continue
    print(f'{n:>6}  {val}')
    print(f'        bits={val.bit_length()}  hex={hex(val)[:40]}...')
    print(f'        mod p = {val % P}')
    for k, c in known.items():
        if val == c: print(f'        *** EQUALS {k}')
        if val % P == c % P and val != c: print(f'        *** == {k} (mod p)')
    if val.bit_length() > 256:
        q, r = divmod(val, P)
        print(f'        = {q} * p + {r}')
        print(f'        quotient bits {q.bit_length()}')
    # small factors
    f, m = [], val
    for d in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97):
        while m % d == 0: f.append(d); m //= d
    print(f'        small factors {f if f else "none"}; cofactor bits {m.bit_length()}')
