"""Map the instance's curve to standard secp256k1 and probe the discrete log."""
import json, time, sys
import dlib as L
P = L.P
B = 64019533680030876408443198762210829058751700634554282185987325820393598524794
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
inv = lambda z: pow(z % P, P - 2, P)

d = json.load(open('ecdlp.json'))
shift = int(d['shift'])
order = d['order']
T = (int(d['target'][0]), int(d['target'][1]))


def onc(pt, b=B):
    x, y = pt
    return (y * y - x ** 3 - b) % P == 0


def add(p1, p2, b=B):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        lam = 3 * x1 * x1 * inv(2 * y1) % P
    else:
        lam = (y2 - y1) * inv(x2 - x1) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(k, pt, b=B):
    r = None
    a = pt
    while k:
        if k & 1:
            r = add(r, a, b)
        a = add(a, a, b)
        k >>= 1
    return r


# --- ladder root ---
t = json.load(open('table.json'))
ent = {int(bb): [C % P for a_, x, C in t[str(bb)]] for bb in map(int, t)}
def toshort(x, y):
    return ((x + shift) % P, y)
pts = {}
for bb, cs in ent.items():
    for x, y in ((cs[0], cs[1]), (cs[1], cs[0])):
        X, Y = toshort(x, y)
        if onc((X, Y)):
            pts[bb] = (X, Y)
P0 = pts[order[0]]
print('ladder root bit x_%d  P0 = %s' % (order[0], P0))
print('P0 on curve:', onc(P0), '  n*P0 = O:', mul(N, P0) is None)
print('target T =', T, ' on curve:', onc(T))

# --- isomorphism to y^2 = x^3 + 7 :  x -> x/u^2, y -> y/u^3 with u^6 = B/7 ---
z = B * inv(7) % P
s = pow(z, (P + 1) // 4, P)
assert s * s % P == z
m = (P - 1) // 3
u = None
for root in (s, (P - s) % P):
    for k in range(9):
        ee = 2 * m + 1 + k * (P - 1)
        if ee % 3 == 0:
            c = pow(root, ee // 3, P)
            if pow(c, 6, P) == z:
                u = c
                break
    if u:
        break
print('u found:', u is not None, ' u^6 == B/7 :', u is not None and pow(u, 6, P) == z)
if u is None:
    u = 1
u2, u3 = u * u % P, pow(u, 3, P)
def tosecp(pt):
    x, y = pt
    return (x * inv(u2) % P, y * inv(u3) % P)
G0 = tosecp(P0)
T0 = tosecp(T)
print('P0 -> secp:', G0, ' on y^2=x^3+7:', (G0[1] ** 2 - G0[0] ** 3 - 7) % P == 0)
print('T  -> secp:', T0, ' on y^2=x^3+7:', (T0[1] ** 2 - T0[0] ** 3 - 7) % P == 0)
print('P0 is the standard secp256k1 generator G:', G0 == (GX, GY))
if G0 != (GX, GY):
    # maybe with the other 6th root branch (y sign / u * zeta)
    print('  -x variant:', (G0[0], (P - G0[1]) % P) == (GX, GY))

# --- BSGS for a small discrete log  T = k*P0, k < M^2 ---
LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 20
print(f'BSGS with m = {LIM} (covers k < {LIM*LIM} ~ 2^{2*LIM.bit_length()-2})', flush=True)
t0 = time.time()
baby = {}
cur = None
for j in range(LIM):
    key = cur[0] if cur else -1
    baby.setdefault(key, j)
    cur = add(cur, P0)
print(f'  baby steps done {time.time()-t0:.0f}s', flush=True)
mP = mul(LIM, P0)
negmP = (mP[0], (P - mP[1]) % P)
g = T
found = None
for i in range(LIM + 1):
    key = g[0] if g else -1
    if key in baby:
        j = baby[key]
        k = (i * LIM + j) % N
        if mul(k, P0) == T:
            found = k
            break
    g = add(g, negmP)
    if i % 50000 == 0:
        print(f'   giant {i} {time.time()-t0:.0f}s', flush=True)
print('discrete log found:', found)
if found is not None:
    json.dump({'k': str(found)}, open('dlog_k.json', 'w'))
