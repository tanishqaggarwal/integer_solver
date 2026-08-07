"""Drive A = x_35389 and B = x_6671 to 0 mod p using x3 = x_22162 and y3 = x_30213."""
import json, sys, time
import dlib as L
import engine2 as E
import adv3
P = L.P
A_, B_ = 35389, 6671
X3, Y3 = 22162, 30213

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json'))
print('base', st0.score, st0.nz())


def AB(st):
    return st.v[A_] % P, st.v[B_] % P


def probe(st0, u, d):
    st = st0.clone()
    st.apply({u: st.v[u] + d})
    adv3.sweep(st, rounds=6)
    return AB(st), st.score


A0, B0 = AB(st0)
print('A0', A0)
print('B0', B0)
for u in (X3, Y3):
    (a1, b1), s1 = probe(st0, u, 1)
    (a2, b2), s2 = probe(st0, u, 2)
    da1, db1 = (a1 - A0) % P, (b1 - B0) % P
    da2, db2 = (a2 - A0) % P, (b2 - B0) % P
    lin = (da2 == 2 * da1 % P) and (db2 == 2 * db1 % P)
    print(f'x_{u}: dA={da1} dB={db1} linear={lin} score(+1)={s1}')

# solve the 2x2 system over GF(p):  A0 + dA_x3*t + dA_y3*s = 0 ; B0 + ... = 0
(a1, b1), _ = probe(st0, X3, 1)
(a2, b2), _ = probe(st0, Y3, 1)
m11, m21 = (a1 - A0) % P, (b1 - B0) % P     # d/dx3
m12, m22 = (a2 - A0) % P, (b2 - B0) % P     # d/dy3
det = (m11 * m22 - m12 * m21) % P
print('det', det)
if det == 0:
    print('singular; falling back to sequential solve')
    sys.exit(0)
di = pow(det, P - 2, P)
rhs1, rhs2 = (-A0) % P, (-B0) % P
t = (m22 * rhs1 - m12 * rhs2) % P * di % P
s = (-m21 * rhs1 + m11 * rhs2) % P * di % P
print('delta x3 =', t)
print('delta y3 =', s)

st = st0.clone()
st.apply({X3: st.v[X3] + t, Y3: st.v[Y3] + s})
adv3.sweep(st, rounds=8)
a, b = AB(st)
print('after joint move: A =', a, ' B =', b, ' score =', st.score, ' nz =', st.nz())
if st.score > 39013:
    L.save(st.v, f'D_ec_{st.score}.json')
    print('saved D_ec_%d.json' % st.score)
L.save(st.v, 'D_ec_last.json')
