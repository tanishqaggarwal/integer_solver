"""Finish the (x_22152, x_33462) construction: after fixing x1,y1 mod p, use the
remaining freedom K = target + p*u to also meet 6672769 | x_25739, then close handles."""
import json, sys, time
import dlib as L
import engine2 as E
import adv3
import rad
P = L.P
K1, K2 = 22152, 33462
A_, B_ = 35389, 6671
TARG = [11150, 25739, 37758]
M2 = 6672769

st0 = E.St(L.load('D_adv.json'))
NX1 = int(sys.argv[1]) if len(sys.argv) > 1 else None
NY1 = int(sys.argv[2]) if len(sys.argv) > 2 else None
if NX1 is None:
    NX1 = 20302955751113177691132960011219991444785130617995423281601414462835238472546
    NY1 = 4531249068709477613185164105669741036354237152756954144434674493737552368539


def run(u, w, sweeps=12):
    st = st0.clone()
    st.apply({K1: NX1 + P * u, K2: NY1 + P * w})
    adv3.sweep(st, rounds=sweeps)
    return st


def tv(st):
    return [st.v[x] for x in TARG]


v00 = tv(run(0, 0))
print('at (0,0):', [x % P for x in v00], '   x_25739 % M2 =', v00[1] % M2)
pts = {}
for (u, w) in [(1, 0), (2, 0), (0, 1), (0, 2), (1, 1), (3, 0)]:
    pts[(u, w)] = tv(run(u, w))
    print(f'  ({u},{w}) x_25739 % M2 = {pts[(u,w)][1] % M2}   x_11150%p={pts[(u,w)][0]%P} x_37758%p={pts[(u,w)][2]%P}')

# fit x_25739 mod M2 as polynomial in u (w=0)
f = [v00[1] % M2, pts[(1, 0)][1] % M2, pts[(2, 0)][1] % M2, pts[(3, 0)][1] % M2]
d1 = [(f[i + 1] - f[i]) % M2 for i in range(3)]
d2 = [(d1[i + 1] - d1[i]) % M2 for i in range(2)]
print('differences in u:', d1, d2)
sol = None
if d2[0] % M2 == 0 and d2[1] % M2 == 0:
    # affine: f0 + d*u == 0 mod M2
    from math import gcd
    d = d1[0] % M2
    g = gcd(d, M2)
    if f[0] % g == 0:
        m = M2 // g
        u = (-(f[0] // g)) * pow((d // g) % m, -1, m) % m
        sol = (u, 0)
print('u solution:', sol)
if sol is None:
    # quadratic in u: brute force over M2 is too big; try w instead
    g = [v00[1] % M2, pts[(0, 1)][1] % M2, pts[(0, 2)][1] % M2]
    e1 = [(g[i + 1] - g[i]) % M2 for i in range(2)]
    print('differences in w:', e1)
    from math import gcd
    d = e1[0] % M2
    gg = gcd(d, M2)
    if g[0] % gg == 0:
        m = M2 // gg
        w = (-(g[0] // gg)) * pow((d // gg) % m, -1, m) % m
        sol = (0, w)
    print('w solution:', sol)
if sol is None:
    sys.exit(1)
u, w = sol
st = run(u, w, sweeps=14)
print('after CRT move: score', st.score, 'nz', st.nz())
print('  x_25739 % (M2*P) =', st.v[25739] % (M2 * P))


def fixatom(st, c):
    c0 = st.av[c]
    if c0 == 0:
        return True
    for uu in sorted(rad.free_knobs(c, st.v)):
        b = st.v[uu]
        r = st.apply({uu: b + 1})
        s = st.av[c] - c0
        st.revert(r)
        if s == 0 or c0 % s:
            continue
        r = st.apply({uu: b - c0 // s})
        if st.av[c] == 0:
            print(f'   fixed a{c} via x_{uu} -> {st.score}')
            return True
        st.revert(r)
    return False


for it in range(6):
    nz = st.nz()
    prog = False
    for c in [19297, 19299, 30984, 36185, 40812]:
        if st.av[c] != 0 and fixatom(st, c):
            prog = True
    if not prog:
        break
print('FINAL', st.score, st.nz())
L.save(st.v, f'D_ec4_{st.score}.json')
print('saved D_ec4_%d.json' % st.score)
