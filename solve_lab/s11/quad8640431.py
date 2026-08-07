"""At the 39,018 state only a26719/a26721/a26723 remain.  They need
      8640431 * P  |  x_12000        (a26721,a26723 need only P | .)
   Shifting x_31339 by k*P keeps x_3719, x_25118 == 0 mod P and moves
      gamma(k) = x_12000/P  quadratically in k.  Solve gamma(k) == 0 mod 8640431.
"""
import sys, os, json, time, random, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
M = 8640431
Q = [53, 163027]
assert Q[0] * Q[1] == M


def sqrt_mod(a, q):
    a %= q
    if a == 0:
        return [0]
    if pow(a, (q - 1) // 2, q) != 1:
        return []
    if q % 4 == 3:
        r = pow(a, (q + 1) // 4, q)
        return sorted({r, q - r})
    # Tonelli-Shanks
    s, e = q - 1, 0
    while s % 2 == 0:
        s //= 2
        e += 1
    n = 2
    while pow(n, (q - 1) // 2, q) != q - 1:
        n += 1
    x = pow(a, (s + 1) // 2, q)
    b = pow(a, s, q)
    g = pow(n, s, q)
    r = e
    while True:
        t, m = b, 0
        for m in range(r):
            if t == 1:
                break
            t = t * t % q
        if m == 0:
            return sorted({x, q - x})
        gs = pow(g, 1 << (r - m - 1), q)
        g = gs * gs % q
        x = x * gs % q
        b = b * g % q
        r = m


def quad_roots(a, b, c, q):
    """roots of a k^2 + b k + c == 0 mod q (q odd prime)"""
    a %= q
    b %= q
    c %= q
    if a == 0:
        if b == 0:
            return list(range(q)) if c == 0 else []
        return [(-c) * pow(b, -1, q) % q]
    disc = (b * b - 4 * a * c) % q
    out = []
    for r in sqrt_mod(disc, q):
        out.append((-b + r) * pow(2 * a, -1, q) % q)
    return sorted(set(out))


def load(name):
    v = [0] * L.NVARS
    d = json.load(open(os.path.join(HERE, 'data', name)))
    for k, x in d.items():
        v[int(k)] = int(x)
    fw.forward(v)
    return v


v = load('three.json')
print("state bad:", fw.bad_checks(v), " failing:", len(L.failing_eqs(L.all_atom_values(v))))
print(f"  x3719%P==0:{v[3719]%P==0}  x25118%P==0:{v[25118]%P==0}")

# gamma(k) by interpolation at k = 0,1,2  (quadratic)
base31339 = v[31339]
vals = []
for k in range(3):
    v[31339] = base31339 + k * P
    fw.forward(v)
    assert v[3719] % P == 0 and v[25118] % P == 0, k
    vals.append((v[12000] // P) % M)
v[31339] = base31339
fw.forward(v)
g0, g1, g2 = vals
# gamma(k) = c + b k + a k^2
c = g0 % M
a = ((g2 - 2 * g1 + g0) * pow(2, -1, M)) % M
b = (g1 - g0 - a) % M
print(f"  gamma(k) = {a} k^2 + {b} k + {c}  (mod {M});  check k=2: {(a*4+b*2+c)%M == g2}")

sols = []
for r53 in quad_roots(a, b, c, Q[0]):
    for r163 in quad_roots(a, b, c, Q[1]):
        k = (r53 * Q[1] * pow(Q[1], -1, Q[0]) + r163 * Q[0] * pow(Q[0], -1, Q[1])) % M
        sols.append(k)
print(f"  roots mod 53: {quad_roots(a,b,c,Q[0])}  mod 163027: {len(quad_roots(a,b,c,Q[1]))} found")
print(f"  CRT solutions: {len(sols)}")

for k in sols[:6]:
    v[31339] = base31339 + k * P
    fw.forward(v)
    gam = (v[12000] // P) % M
    bad = fw.bad_checks(v)
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"  k={k}: gamma%M={gam} bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:10]}", flush=True)
    if gam == 0:
        for a2 in list(bad):
            for t, d in ([(u, None) for u in L.avars[a2] if L.definer.get(u) is None] +
                         [(t, d) for t, d in deep.handles(v, a2, locked=set())[0]]):
                old = v[t]
                if d is None:
                    x = fw.solve_lin(a2, t, v)
                    if x is None or x == old:
                        continue
                else:
                    bs = fw.evalpoly(L.polys[a2], v)
                    if not d or bs % d:
                        continue
                    x = old - bs // d
                v[t] = x
                fw.forward(v)
                if fw.evalpoly(L.polys[a2], v) == 0:
                    break
                v[t] = old
                fw.forward(v)
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        print(f"     after close: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:10]}", flush=True)
        json.dump({str(i): v[i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', f'quad_k{k}.json'), 'w'))
v[31339] = base31339
fw.forward(v)
