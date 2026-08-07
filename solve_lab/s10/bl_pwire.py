"""bl_pwire: what pins the p-wires that the seven residual checks need?

In frame2 the seven checks reduce to congruences on the FREE p-wires
  x_9118, x_8731, x_1329, x_10903, x_17325, x_9413, x_6418, x_31861
(plus the detached x_642, x_7068, x_28730, x_29854, x_31864).
Find every atom / equation that constrains those wires, and ask which of those
constraints is gated by a boolean (i.e. would vanish if some boolean flipped).
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, BOOLATOM, CANON, F2, pot, FORBID
P = 2**256-2**32-977

w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
av = L.all_atom_values(w)
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
print('residual atom values (as integers, and mod p):')
for a in SEVEN:
    print(f'  a{a}: {str(av[a])[:40]}...  bits={av[a].bit_length()}  mod p = {av[a]%P}')

WIRES = [9118, 8731, 1329, 10903, 17325, 9413, 6418, 31861, 642, 7068, 28730, 29854, 31864, 2099]
print('\nwire values mod p:')
for u in WIRES:
    q, r = divmod(w[u], P)
    print(f'  x_{u:<6} bits={w[u].bit_length():<5} mod p = {r}   (q bits {q.bit_length()})')

print('\nrequired congruences for the seven to vanish (branch x_7075=1):')
print(f'  need P | 5113045*x_9118  -> x_9118 mod p = {w[9118]%P} (need 0)')
print(f'  need P | x_8731          -> x_8731 mod p = {w[8731]%P} (need 0)')
print(f'  a22229: x_7068 - x_2099 - 7376877*x_642, with x_642 = P*x_17325')
print(f'          (x_7068-x_2099) mod p = {(w[7068]-w[2099])%P}')

# --- who else uses these wires? ---
print('\n--- atoms (other than the seven) that mention each wire ---')
for u in (9118, 8731, 1329, 10903, 17325, 9413, 6418, 31861):
    As = [a for a in L.var_atoms[u] if a not in SEVEN]
    eqs = set()
    for a in As: eqs |= set(L.atom2eq.get(a, ()))
    print(f'  x_{u:<6}: {len(As)} atoms, {len(eqs)} equations   atoms={As[:14]}')
    for a in As[:6]:
        print(f'        a{a}: {L.atom_src[a][:110]}')

# --- boolean-gated dependence: for each atom mentioning a wire, does a boolean
#     multiply the wire (so a boolean=0 kills the term)?
print('\n--- boolean gates sitting on those wires ---')
hits = collections.defaultdict(set)
for u in (9118, 8731, 1329, 10903, 17325, 9413, 6418, 31861):
    for a in L.var_atoms[u]:
        for m, c in L.polys[a].items():
            if u in m:
                for z in m:
                    if z != u and z in BOOL:
                        hits[u].add((a, z))
for u, s in hits.items():
    print(f'  x_{u}: gated in {len(s)} monomials by booleans {sorted(set(z for _, z in s))}')
    for a, z in sorted(s)[:8]:
        print(f'      a{a} boolean x_{z} (free={z in CANON.FREE}, val={w[z]})  {L.atom_src[a][:90]}')
if not hits: print('  NONE')
