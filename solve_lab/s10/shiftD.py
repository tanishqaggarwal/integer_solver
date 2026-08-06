"""S10 step 16: can D = x_7068 - x_2099 be shifted by multiples of p?

The residual constraint is  A1 + 7376877*A7 == D  (mod 7376877*p).
D mod p is the hard core congruence, but D mod 7376877 may be movable: shifting
x_7068 by k*p leaves every p-quantised link intact while changing D mod 7376877
(gcd(p,7376877)=1, so k*p covers all residues).  Test it.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
base = L.load(BEST)
bav = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if bav[a])
print('base nonzero atoms:', sorted(BASE_NZ))

for tag, delta in (('+1', 1), ('+p', P), ('+2p', 2 * P), ('+1000p', 1000 * P),
                   ('-3228258p', -3228258 * P)):
    for var in (7068, 2099):
        v = list(base)
        ch, _ = L.ripple(v, {var: base[var] + delta}, block=BLOCK)
        av = L.all_atom_values(v)
        nz = set(a for a in range(L.NA) if av[a])
        extra = sorted(nz - BASE_NZ)
        gone = sorted(BASE_NZ - nz)
        D = v[7068] - v[2099]
        print(f'x_{var} {tag:<12} changed={len(ch):<4} extra_nz={extra} '
              f'gone={gone} D%7376877={D % 7376877} D%p==base:{D % P == (base[7068]-base[2099]) % P}')
