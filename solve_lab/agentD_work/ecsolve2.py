"""Exact joint solve.

At D_adv the only nonzero atoms are
  a19297 = sel*x_11150 + x_4007      (x_4007  = p*h1)          -> need p | x_11150
  a19299 = sel*x_25739 - 6672769*x_29804  (x_29804 = p*h2)     -> need 6672769*p | x_25739
  a30984 = 537773*(sel*x_37758) - x_35605 (x_35605 = p*h3)     -> need p | x_37758
plus two 1-equation shadows.

x_11150, x_25739, x_37758 are exactly affine over Z in t = delta(x_22162)=x3 and
s = delta(x_30213)=y3.  Measure the affine map by probing, then solve the
congruences by CRT and verify by exact evaluation.
"""
import json, sys, time
import dlib as L
import engine2 as E
import adv3
P = L.P
X3, Y3 = 22162, 30213
TARG = [11150, 25739, 37758]

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json'))
print('base', st0.score, st0.nz())
base_t, base_s = st0.v[X3], st0.v[Y3]


def run(t, s, sweeps=8):
    st = st0.clone()
    st.apply({X3: base_t + t, Y3: base_s + s})
    adv3.sweep(st, rounds=sweeps)
    return st


def vals(st):
    return [st.v[u] for u in TARG]


v00 = vals(st0)
v10 = vals(run(1, 0))
v01 = vals(run(0, 1))
a_t = [v10[i] - v00[i] for i in range(3)]
a_s = [v01[i] - v00[i] for i in range(3)]
v11 = vals(run(1, 1))
ok = all(v11[i] == v00[i] + a_t[i] + a_s[i] for i in range(3))
print('affine check (1,1):', ok)
v22 = vals(run(2, 3))
ok2 = all(v22[i] == v00[i] + 2 * a_t[i] + 3 * a_s[i] for i in range(3))
print('affine check (2,3):', ok2)
for i, u in enumerate(TARG):
    print(f'  x_{u}: c0 bits={v00[i].bit_length()} dt bits={a_t[i].bit_length()} ds bits={a_s[i].bit_length()}')

M2 = 6672769

# --- solve mod p:  v00[i] + a_t[i]*t + a_s[i]*s == 0 (mod p) for i = 0,2 ---
def solve2(rows, mod):
    """rows = [(c, ct, cs)] ; solve ct*t + cs*s == -c (mod mod). Returns (t,s) or None."""
    import itertools
    (c1, t1, s1), (c2, t2, s2) = rows
    det = (t1 * s2 - t2 * s1) % mod
    from math import gcd
    if gcd(det, mod) != 1:
        return None
    di = pow(det, -1, mod)
    r1, r2 = (-c1) % mod, (-c2) % mod
    t = (s2 * r1 - s1 * r2) % mod * di % mod
    s = (-t2 * r1 + t1 * r2) % mod * di % mod
    return t, s


rp = solve2([(v00[0], a_t[0], a_s[0]), (v00[2], a_t[2], a_s[2])], P)
print('mod p solution:', rp is not None)
tp, sp = rp
# check the middle row also vanishes mod p
chk = (v00[1] + a_t[1] * tp + a_s[1] * sp) % P
print('x_25739 mod p at that solution:', chk)

# --- now fix t = tp + P*u, s = sp + P*w and solve x_25739 == 0 (mod 6672769) ---
c = v00[1] + a_t[1] * tp + a_s[1] * sp
ct = a_t[1] * P
cs = a_s[1] * P
from math import gcd
g = gcd(gcd(ct % M2, cs % M2), M2)
print('M2 =', M2, ' gcd(ct,cs,M2) =', g, ' c mod g =', c % g)
sol = None
if c % g == 0:
    m = M2 // g
    ctg, csg, cg = (ct // g) % m if g else ct, (cs // g) % m, (c // g) % m
    # solve ctg*u + csg*w == -cg (mod m)
    if gcd(ctg, m) == 1:
        u = (-cg) * pow(ctg, -1, m) % m
        w = 0
        sol = (u, w)
    elif gcd(csg, m) == 1:
        w = (-cg) * pow(csg, -1, m) % m
        u = 0
        sol = (u, w)
    else:
        for w in range(0, min(m, 200000)):
            rr = (-cg - csg * w) % m
            gg = gcd(ctg, m)
            if rr % gg == 0:
                mm = m // gg
                u = (rr // gg) * pow((ctg // gg) % mm, -1, mm) % mm
                sol = (u, w)
                break
print('CRT solution:', sol)
if sol is None:
    sys.exit(1)
u, w = sol
t = tp + P * u
s = sp + P * w
print('t bits', t.bit_length(), 's bits', s.bit_length())
st = run(t, s, sweeps=10)
print('score after move', st.score, 'nz', st.nz())
for i, uu in enumerate(TARG):
    val = st.v[uu]
    print(f'  x_{uu} % p = {val % P}   % (6672769*P) = {val % (M2*P)}')
L.save(st.v, 'D_ec2.json')
print('saved D_ec2.json')
