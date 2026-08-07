"""S10 step 88: read the DELIVERED witness in the cone's coordinates."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
W = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(W)
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
nz = [a for a in range(L.NA) if av[a]]
print(f'delivered witness: failing {len(fail)}  nonzero atoms {nz}')
for a in nz:
    print(f'  a{a:<6} eqs={len(L.atom2eq[a]):>3} val={str(av[a])[:40]}  '
          f'val mod p = {str(av[a] % P)[:40]}')
print()
for t in [7068, 2099, 28730, 642, 9413, 17325, 6418, 9118, 31861, 2081, 4287,
          17499, 28599, 26064]:
    print(f'  x_{t:<6} = {str(v[t])[:46]:<48} (bits {v[t].bit_length()})')
D = v[7068] - v[2099]
print(f'\nD = x_7068 - x_2099 = {str(D)[:50]}  (bits {D.bit_length()})')
print(f'D mod p  = {D % P}')
print(f'D mod 7376877 = {D % 7376877}')
print(f'x_28730 mod p = {v[28730] % P}')
print(f'x_26064 = p ? {v[26064] == P};  x_17499 = p ? {v[17499] == P}')
