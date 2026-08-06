"""S10 step 11: is x_9118 / x_8731 free, or pinned mod p?

clean.py showed that setting them to 0 breaks atoms 7930 and 41512.  Determine
whether those atoms constrain the residue mod p (in which case congruences 2,3
are NOT free) or the exact integer.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')

for var in (9118, 8731):
    print(f'\n########## x_{var} ##########')
    for name, delta in (('+1', 1), ('+p', P), ('+2p', 2 * P), ('+7376877*p', 7376877 * P)):
        v = L.load(BEST)
        base = v[var]
        ch, _ = L.ripple(v, {var: base + delta})
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        fail = L.failing_eqs(av)
        print(f'  delta {name:<12} changed={len(ch):<4} nz_atoms={len(nz):<3} '
              f'failing={len(fail):<4} score={L.NEQ-len(fail)}')
        print(f'      nz = {nz}')
