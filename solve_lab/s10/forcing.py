"""S11 step 54: verify the FORCING CHAIN and the constant mismatch.

Claim: a7930 forces x_24548 == x_12553 ; a3578 pins x_12553 == C3 ; a21617 forces
x_14623 == x_24548 ; the K1 web plus a31672 pins x_14623 == K1 -- and C3 != K1
mod p.  If the chain holds with the selectors on, the instance has NO full
solution in that branch.
"""
import os, sys
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
def pr(a, n=200):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
for a in (7930, 3578, 21617, 31672, 3576):
    print(f'a{a} ({len(L.atom2eq[a])} eqs): {pr(a)}')
# extract the pinned constants
def const_of(a, sel, var):
    """atom = sel*(var - C) - ... : recover C from the coefficient of the sel monomial"""
    c = L.polys[a].get((sel,), 0)
    return -c
C3 = const_of(3578, 2081, 12553)
K1 = const_of(31672, 24601, 33462)
C4 = const_of(3576, 2081, 6418)
print(f'\nC3 = {C3}\n   C3 mod p = {C3 % P}')
print(f'K1 = {K1}\n   K1 mod p = {K1 % P}')
print(f'C4 = {C4}\n   C4 mod p = {C4 % P}')
print(f'\n*** C3 == K1 (mod p) ? {C3 % P == K1 % P}')
print(f'    C3 - K1 mod p = {(C3 - K1) % P}')
print(f'    gcd-style: is C3 == K1 exactly? {C3 == K1}')
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
print(f'\nat the delivered witness:')
for u in (2081, 24601, 12553, 24548, 14623, 33462, 6418):
    print(f'  x_{u:<6} = {str(v[u])[:40]}   mod p = {str(v[u] % P)[:40]}')
print(f'\n  x_12553 == C3 (mod p)? {v[12553] % P == C3 % P}')
print(f'  x_33462 == K1 (mod p)? {v[33462] % P == K1 % P}')
print(f'  x_24548 == x_12553 (mod p)? {v[24548] % P == v[12553] % P}')
print(f'  x_14623 == x_24548 (mod p)? {v[14623] % P == v[24548] % P}')
print(f'  x_14623 == x_33462 (mod p)? {v[14623] % P == v[33462] % P}')
