"""Are the region's four multiplier constants equal to p?"""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
v = G.V0
for u in [17499, 22665, 28961, 28599, 21279, 7075]:
    x = v[u]
    print(f'x_{u} = {abs(x).bit_length()}b   == p: {x == P}   == -p: {x == -P}   '
          f'p | x: {x % P == 0 if x else "0"}   x = {str(x)[:40]}')
print('p =', str(P)[:40], f'({P.bit_length()}b)')
# the four blocked numerators mod p
print()
R = G.R0 + [23618]
Pv = G.private_vars(R)
const, cols = G.build_model(R, Pv, v)
Eqs, rows = G.eq_system(R, Pv, const, cols)
print('columns of each private knob, divisibility by p:')
for u in Pv:
    cs = [c for c in cols[u].values()]
    print(f'  x_{u}: entries {[abs(c).bit_length() for c in cs]}  all divisible by p: '
          f'{all(c % P == 0 for c in cs)}')
